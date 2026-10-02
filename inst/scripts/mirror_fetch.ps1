<#
Copy network files or directories into the EBVhelpR local data mirror and
record where they came from. Windows counterpart of EBVhelpR::ebv_mirror_fetch():
same mirror layout, same MANIFEST.tsv. Use it for big files, where robocopy over
SMB (~110 MB/s measured) beats R reading through a WSL 9p mount (~3 MB/s).

  powershell -File mirror_fetch.ps1 <source> [<source> ...]
  powershell -File mirror_fetch.ps1 -NoHash <big source>

A source is a mapped-drive or UNC path to a file or a directory (copied
recursively). Mapped drives are resolved to UNC, so a file lands in the same
place whichever letter reached it (the lab share is L: on one machine, Z: on
another). The mirror path is cache\<UNC without the leading \\>\..., which is the
`mirror_subdir` of the matching location in EBVhelpR's data_locations.csv:

  L:\Centers\CBSR\PI\Volaric\x.csv  ->  <mirror>\cache\files.med.uvm.edu\shared\Centers\CBSR\PI\Volaric\x.csv

Mirror root: $env:EBVHELPER_MIRROR_DIR_WIN, else %USERPROFILE%\EBV_data_mirror.

COPY ONLY. Nothing at the source is touched and an existing mirror file is never
overwritten: if it differs from the source it is left alone and a row is added
to <mirror>\cleanup_targets.tsv. Each copied file appends one manifest row
(source, local, bytes, source mtime, md5 unless -NoHash, fetch time, machine).
#>
param(
  [Parameter(Mandatory = $true, ValueFromRemainingArguments = $true)][string[]]$Source,
  [switch]$NoHash
)
$ErrorActionPreference = 'Stop'
$mirror = if ($env:EBVHELPER_MIRROR_DIR_WIN) { $env:EBVHELPER_MIRROR_DIR_WIN }
          elseif ($env:P1_MIRROR_WIN) { $env:P1_MIRROR_WIN }
          else { Join-Path $env:USERPROFILE 'EBV_data_mirror' }
$manifest = Join-Path $mirror 'MANIFEST.tsv'
$cleanup  = Join-Path $mirror 'cleanup_targets.tsv'
if (-not (Test-Path $manifest)) {
  New-Item -ItemType Directory -Force $mirror | Out-Null
  "source`tlocal`tbytes`tsource_mtime`tmd5`tfetched`tmachine" | Out-File -Encoding utf8 $manifest
}
function Stamp([datetime]$t) { $t.ToUniversalTime().ToString('yyyy-MM-ddTHH:mm:ssZ') }

function To-Unc([string]$p) {
  $full = [IO.Path]::GetFullPath($p)
  if ($full -match '^([A-Za-z]):\\(.*)$') {
    $drv = Get-PSDrive -Name $Matches[1] -ErrorAction SilentlyContinue
    if ($drv -and $drv.DisplayRoot) { return ($drv.DisplayRoot.TrimEnd('\') + '\' + $Matches[2]) }
  }
  return $full
}

function Log-Cleanup([string]$path, [string]$reason, [string]$keep) {
  if (-not (Test-Path $cleanup)) { "logged_utc`tmachine`tlocation`tpath`tbytes`treason`tkeep" | Out-File -Encoding utf8 $cleanup }
  if (Select-String -LiteralPath $cleanup -SimpleMatch "`t$path`t" -Quiet) { return }
  "$(Stamp (Get-Date))`t$env:COMPUTERNAME`tmirror`t$path`t$((Get-Item -LiteralPath $path).Length)`t$reason`t$keep" |
    Out-File -Append -Encoding utf8 $cleanup
}

function Fetch-File([IO.FileInfo]$item) {
  $unc = To-Unc $item.FullName
  if ($unc -notmatch '^\\\\') { throw "not a network path (nothing to mirror): $unc" }
  $dest = Join-Path (Join-Path $mirror 'cache') $unc.TrimStart('\')
  if (Test-Path -LiteralPath $dest) {
    $d = Get-Item -LiteralPath $dest
    $same = $d.Length -eq $item.Length -and [math]::Abs(($d.LastWriteTimeUtc - $item.LastWriteTimeUtc).TotalSeconds) -le 2
    if (-not $same) { Log-Cleanup $dest 'mirror copy differs from its source; not overwritten (copy-only policy)' $unc; return 'differs' }
    if (Select-String -LiteralPath $manifest -SimpleMatch "$unc`t" -Quiet) { return 'unchanged' }
    # copied by hand earlier: only the manifest row is missing
  } else {
    $null = robocopy (Split-Path $item.FullName) (Split-Path $dest) $item.Name /J /NP /NJH /NJS /R:2 /W:5 /COPY:DT /XC /XN /XO
    if ($LASTEXITCODE -ge 8) { throw "robocopy failed ($LASTEXITCODE) for $($item.FullName)" }
  }
  $md5 = if ($NoHash) { '' } else { (Get-FileHash -Algorithm MD5 -LiteralPath $dest).Hash.ToLower() }
  "$unc`t$dest`t$($item.Length)`t$(Stamp $item.LastWriteTime)`t$md5`t$(Stamp (Get-Date))`t$env:COMPUTERNAME" |
    Out-File -Append -Encoding utf8 $manifest
  return 'copied'
}

foreach ($s in $Source) {
  if (Test-Path -LiteralPath $s -PathType Container) {
    $files = Get-ChildItem -LiteralPath $s -Recurse -File -Force | Where-Object { $_.Name -ne '.ebv_fetch_complete' }
    $res = $files | ForEach-Object { Fetch-File $_ }
    $n = ($res | Group-Object | ForEach-Object { "$($_.Count) $($_.Name)" }) -join ', '
    if (-not ($res -contains 'differs')) {
      $destDir = Join-Path (Join-Path $mirror 'cache') (To-Unc (Get-Item -LiteralPath $s).FullName).TrimStart('\')
      "source: $s`r`nfiles: $($files.Count)`r`ncompleted: $(Stamp (Get-Date))" | Out-File -Encoding ascii (Join-Path $destDir '.ebv_fetch_complete')
    }
    Write-Output "$s : $n"
  } elseif (Test-Path -LiteralPath $s -PathType Leaf) {
    Write-Output "$(Fetch-File (Get-Item -LiteralPath $s))  $s"
  } else { throw "not reachable: $s" }
}
