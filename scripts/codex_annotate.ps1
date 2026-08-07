#!/usr/bin/env pwsh
# CellTypePilot — Quick annotation script for Codex CLI (Windows PowerShell)
# Usage: .\scripts\codex_annotate.ps1 -Input <path> -ClusterKey <key> [-OutputDir <dir>] [-Tissue <tissue>]
#
# This script wraps the full CellTypePilot annotation pipeline into a single
# command that Codex can execute without multi-step orchestration.

param(
    [Parameter(Mandatory=$true)]
    [string]$Input,

    [Parameter(Mandatory=$true)]
    [string]$ClusterKey,

    [string]$OutputDir = ".\ctp_output",

    [string]$Tissue = ""
)

Write-Host "=== CellTypePilot — Codex Quick Annotate ===" -ForegroundColor Cyan
Write-Host ""

# Step 1: Environment check
Write-Host "[1/4] Checking environment..." -ForegroundColor Blue
try {
    $null = & celltypepilot doctor 2>&1
} catch {
    Write-Host "ERROR: celltypepilot not found. Run: pip install -e ." -ForegroundColor Red
    exit 1
}

# Step 2: Inspect data
Write-Host "[2/4] Inspecting data..." -ForegroundColor Blue
$inspectJson = & celltypepilot inspect --input $Input --cluster-key $ClusterKey --json 2>$null | ConvertFrom-Json
$species = if ($inspectJson.species) { $inspectJson.species } else { "human" }

if (-not $Tissue) {
    $Tissue = if ($inspectJson.tissue) { $inspectJson.tissue } else { "general" }
}

$embedding = ""
if ($inspectJson.embedding_keys) {
    foreach ($k in $inspectJson.embedding_keys) {
        if ($k -match "umap") {
            $embedding = $k
            break
        }
    }
    if (-not $embedding -and $inspectJson.embedding_keys.Count -gt 0) {
        $embedding = $inspectJson.embedding_keys[0]
    }
}

Write-Host "  Species: $species"
Write-Host "  Tissue:  $Tissue"
Write-Host "  Embedding: $(if ($embedding) { $embedding } else { 'auto' })"
Write-Host ""

# Step 3: Run annotation
Write-Host "[3/4] Running annotation pipeline..." -ForegroundColor Blue
$args = @("annotate", "--input", $Input, "--cluster-key", $ClusterKey, "--output", $OutputDir, "--species", $species, "--tissue", $Tissue)
if ($embedding) {
    $args += @("--embedding-key", $embedding)
}

& celltypepilot @args

# Step 4: Summary
Write-Host ""
Write-Host "[4/4] Done! Output files:" -ForegroundColor Green
Write-Host "  $OutputDir\data.annotated.h5ad"
Write-Host "  $OutputDir\evidence_table.csv"
Write-Host "  $OutputDir\report_draft.html"
Write-Host "  $OutputDir\manifest.json"
Write-Host "  $OutputDir\methodology_draft.txt"
if (Test-Path "$OutputDir\figures") {
    $figCount = (Get-ChildItem "$OutputDir\figures" -File).Count
    Write-Host "  $OutputDir\figures\ ($figCount files)"
}
Write-Host ""
Write-Host "=== Annotation complete ===" -ForegroundColor Green
