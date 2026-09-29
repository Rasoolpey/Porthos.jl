<#
.SYNOPSIS
    Install everything Porthos.jl needs on a Windows machine, in one go.

.DESCRIPTION
    Safe to re-run: each step checks what is already there and skips it.

      1. juliaup (the official Julia installer, via winget), user-level.
      2. Julia 1.12, pinned for this directory with a juliaup override (the global
         default Julia is left alone).
      3. The Julia packages from Project.toml / Manifest.toml (Pkg.instantiate) and
         precompilation.
      4. The parity pack (lazy artifact in Artifacts.toml). This fails harmlessly until the
         GitHub release parity-pack-v1 exists; the step then only warns.
      5. The Python venv for the parity-pack generator (parity/generate/.venv) with the
         pinned parity/generate/requirements.txt. Only needed to regenerate the pack; skip
         with -SkipPython.
      6. A smoke test: load Porthos.

    Not installed yet: a C++ compiler and SUNDIALS, which PHPS needs for its compiled DAE
    runs (the P7 reference trajectories). See TODO.md.

.PARAMETER SkipPython
    Do not create the parity-pack generator venv.

.PARAMETER RunTests
    Run the test suite at the end (about 15 s once precompiled).

.EXAMPLE
    powershell -ExecutionPolicy Bypass -File scripts\setup.ps1
    powershell -ExecutionPolicy Bypass -File scripts\setup.ps1 -RunTests
#>
[CmdletBinding()]
param(
    [switch]$SkipPython,
    [switch]$RunTests
)

$ErrorActionPreference = 'Stop'
$Root = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$JuliaChannel = '1.12'
$PythonVersion = '3.13'

function Step($msg) { Write-Host ""; Write-Host "==> $msg" -ForegroundColor Cyan }
function Ok($msg) { Write-Host "    $msg" -ForegroundColor Green }
function Warn($msg) { Write-Host "    $msg" -ForegroundColor Yellow }

function Refresh-Path {
    $env:Path = [Environment]::GetEnvironmentVariable('Path', 'User') + ';' +
                [Environment]::GetEnvironmentVariable('Path', 'Machine')
}

# Run a native command for its exit code only, tolerating output on stderr (Windows
# PowerShell 5.1 turns redirected stderr into an error under ErrorActionPreference Stop).
function Test-Native([string]$exe, [string[]]$arguments) {
    $old = $ErrorActionPreference
    $ErrorActionPreference = 'Continue'
    try { & $exe @arguments *> $null; return ($LASTEXITCODE -eq 0) }
    catch { return $false }
    finally { $ErrorActionPreference = $old }
}

# Run a native command; fail the script if it exits non-zero.
function Invoke-Native([string]$exe, [string[]]$arguments) {
    & $exe @arguments
    if ($LASTEXITCODE -ne 0) { throw "$exe $($arguments -join ' ') failed (exit $LASTEXITCODE)" }
}

# 1. juliaup ---------------------------------------------------------------------------
Step 'juliaup'
if (-not (Get-Command juliaup -ErrorAction SilentlyContinue)) {
    if (-not (Get-Command winget -ErrorAction SilentlyContinue)) {
        throw 'winget not found. Install juliaup by hand: https://julialang.org/install/'
    }
    Invoke-Native winget @('install', '--id', '9NJNWW8PVKMN', '-e',
                           '--accept-package-agreements', '--accept-source-agreements',
                           '--disable-interactivity')
    Refresh-Path
    if (-not (Get-Command juliaup -ErrorAction SilentlyContinue)) {
        throw 'juliaup was installed but is not on PATH yet; open a new terminal and re-run.'
    }
    Ok 'installed'
} else {
    Ok 'already installed'
}

# 2. Julia 1.12, pinned for this directory ---------------------------------------------
Step "Julia $JuliaChannel"
$old = $ErrorActionPreference; $ErrorActionPreference = 'Continue'
$installed = (& juliaup status 2>&1 | Out-String)
$ErrorActionPreference = $old
if ($installed -notmatch [regex]::Escape($JuliaChannel)) {
    Invoke-Native juliaup @('add', $JuliaChannel)
} else {
    Ok "channel $JuliaChannel already installed"
}
Invoke-Native juliaup @('override', 'set', '--path', $Root, $JuliaChannel)
Push-Location $Root
try {
    $v = (& julia --version)
    Ok "$v (juliaup override for $Root)"
} finally { Pop-Location }

# 3. Julia packages -----------------------------------------------------------------------
Step 'Julia packages (Pkg.instantiate + precompile; the first run takes several minutes)'
Push-Location $Root
try {
    Invoke-Native julia @('--project=.', '-e', 'using Pkg; Pkg.instantiate(); Pkg.precompile()')
    Ok 'done'
} finally { Pop-Location }

# 4. Parity pack ------------------------------------------------------------------------
Step 'Parity pack artifact'
Push-Location $Root
try {
    # (No double quotes inside the Julia code: PowerShell 5.1 and 7 pass them differently.)
    $ok = Test-Native julia @('--project=.', '-e',
        'using Pkg, Porthos; Pkg.Artifacts.ensure_artifact_installed(Porthos.PARITY_ARTIFACT, Porthos.ARTIFACTS_TOML); load_parity_pack()')
    if ($ok) {
        Ok 'installed and verified'
    } else {
        Warn 'Could not install the parity pack (is the release parity-pack-v1 uploaded?).'
        Warn 'The parity tests need it; see parity/README.md.'
    }
} finally { Pop-Location }

# 5. Python venv for the parity-pack generator ------------------------------------------
if ($SkipPython) {
    Step 'Python venv: skipped (-SkipPython)'
} else {
    Step "Python $PythonVersion venv for the parity-pack generator"
    $venv = Join-Path $Root 'parity\generate\.venv'
    $py = Join-Path $venv 'Scripts\python.exe'
    if (-not (Test-Path $py)) {
        $havePy = [bool](Get-Command py -ErrorAction SilentlyContinue) -and
                  (Test-Native py @("-$PythonVersion", '-c', 'import sys'))
        if (-not $havePy) {
            if (-not (Get-Command winget -ErrorAction SilentlyContinue)) {
                throw "Python $PythonVersion not found and winget is unavailable."
            }
            Invoke-Native winget @('install', '--id', "Python.Python.$PythonVersion", '-e',
                                   '--accept-package-agreements', '--accept-source-agreements',
                                   '--disable-interactivity')
            Refresh-Path
        }
        Invoke-Native py @("-$PythonVersion", '-m', 'venv', $venv)
        Ok "created $venv"
    } else {
        Ok "venv exists: $venv"
    }
    Invoke-Native $py @('-m', 'pip', 'install', '--quiet', '--upgrade', 'pip')
    Invoke-Native $py @('-m', 'pip', 'install', '--quiet', '-r',
                        (Join-Path $Root 'parity\generate\requirements.txt'))
    Ok 'requirements installed'
}

# 6. Smoke test / tests ------------------------------------------------------------------
Push-Location $Root
try {
    if ($RunTests) {
        Step 'Test suite'
        Invoke-Native julia @('--project=.', 'test/runtests.jl')
    } else {
        Step 'Smoke test'
        Invoke-Native julia @('--project=.', '-e',
            'using Porthos; c = load_case(ARGS[1]); display(c); display(solve_powerflow(c))',
            'cases/IEEE39Bus_PF/system_phtrue.json')
    }
} finally { Pop-Location }

Write-Host ""
Write-Host 'Porthos setup complete.' -ForegroundColor Green
