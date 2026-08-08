"""Environment doctor — check dependencies and report capability level."""

from __future__ import annotations

import importlib
import shutil
import sys
from dataclasses import dataclass, field

from rich.console import Console
from rich.table import Table

console = Console()


@dataclass
class DependencyStatus:
    name: str
    installed: bool
    version: str = ""
    required: bool = True
    note: str = ""


@dataclass
class DoctorReport:
    python_version: str = ""
    python_ok: bool = False
    dependencies: list[DependencyStatus] = field(default_factory=list)
    optional_deps: list[DependencyStatus] = field(default_factory=list)
    capabilities: dict[str, str] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    mcp_status: dict[str, bool] = field(default_factory=dict)


def check_python() -> tuple[bool, str]:
    """Check Python version."""
    ver = sys.version_info
    version_str = f"{ver.major}.{ver.minor}.{ver.micro}"
    ok = ver.major == 3 and ver.minor >= 10
    return ok, version_str


def check_dependency(name: str, required: bool = True) -> DependencyStatus:
    """Check if a Python package is installed."""
    try:
        try:
            from importlib.metadata import version as pkg_version

            version = pkg_version(name)
        except Exception:
            mod = importlib.import_module(name)
            version = getattr(mod, "__version__", "unknown")
        return DependencyStatus(name=name, installed=True, version=version, required=required)
    except ImportError:
        return DependencyStatus(
            name=name,
            installed=False,
            required=required,
            note="NOT INSTALLED" if required else "optional",
        )


def run_doctor() -> DoctorReport:
    """Run full environment check."""
    report = DoctorReport()

    # Python version
    report.python_ok, report.python_version = check_python()
    if not report.python_ok:
        report.warnings.append(f"Python {report.python_version} detected; >= 3.10 required.")

    # Core dependencies
    core_deps = [
        "anndata",
        "scanpy",
        "numpy",
        "pandas",
        "scipy",
        "matplotlib",
        "seaborn",
        "typer",
        "rich",
    ]
    for dep in core_deps:
        report.dependencies.append(check_dependency(dep, required=True))

    # Optional dependencies
    opt_deps = [
        ("python_docx", "docx export", False),
        ("cupy", "GPU acceleration", False),
        ("scvi_tools", "reference mapping (Phase 2)", False),
        ("decoupler", "pathway scoring (Phase 2)", False),
        ("flask", "Web Inspector", False),
        ("rpy2", "Seurat .rds support", False),
    ]
    for dep_name, desc, req in opt_deps:
        status = check_dependency(dep_name, required=req)
        status.note = desc
        report.optional_deps.append(status)

    # External tools
    if shutil.which("pixi"):
        report.capabilities["pixi"] = "available"
    else:
        report.capabilities["pixi"] = "not found (optional — for isolated env)"

    # MCP / Literature integration check
    try:
        from .literature import check_mcp_availability

        report.mcp_status = check_mcp_availability()
    except Exception:
        report.mcp_status = {"pubmed_direct": False}

    # Determine capability level
    all_core_installed = all(d.installed for d in report.dependencies)
    if all_core_installed:
        report.capabilities["tier"] = "full"
        report.capabilities["description"] = (
            "All core dependencies met. Core local plugin functionality available."
        )
    else:
        missing = [d.name for d in report.dependencies if not d.installed]
        report.capabilities["tier"] = "degraded"
        report.capabilities["description"] = (
            f"Missing core dependencies: {', '.join(missing)}. Run: pip install celltypepilot"
        )

    return report


def format_doctor_report(report: DoctorReport) -> str:
    """Format doctor report as rich table string."""
    lines = []
    lines.append("=" * 60)
    lines.append("CellTypePilot — Environment Check (doctor)")
    lines.append("=" * 60)
    lines.append("")

    # Python
    status = "[OK]" if report.python_ok else "[FAIL]"
    lines.append(f"Python: {report.python_version} {status}")
    lines.append("")

    # Core deps table
    table = Table(title="Core Dependencies")
    table.add_column("Package", style="cyan")
    table.add_column("Status", style="bold")
    table.add_column("Version", style="dim")

    for dep in report.dependencies:
        status_str = (
            f"[green]OK {dep.version}[/green]" if dep.installed else "[red]X NOT INSTALLED[/red]"
        )
        table.add_row(dep.name, status_str, dep.version if dep.installed else "")

    # Use rich to render
    from io import StringIO

    str_io = StringIO()
    temp_console = Console(file=str_io, force_terminal=False)
    temp_console.print(table)
    lines.append(str_io.getvalue())

    # Optional deps
    if report.optional_deps:
        lines.append("")
        opt_table = Table(title="Optional Dependencies")
        opt_table.add_column("Package", style="cyan")
        opt_table.add_column("Status")
        opt_table.add_column("Note", style="dim")

        for dep in report.optional_deps:
            if dep.installed:
                status_str = f"[green]OK {dep.version}[/green]"
            else:
                status_str = "[yellow]-- not installed[/yellow]"
            opt_table.add_row(dep.name, status_str, dep.note)

        str_io2 = StringIO()
        temp_console2 = Console(file=str_io2, force_terminal=False)
        temp_console2.print(opt_table)
        lines.append(str_io2.getvalue())

    # Capabilities
    lines.append("")
    lines.append(f"Capability tier: {report.capabilities.get('tier', 'unknown')}")
    lines.append(f"  {report.capabilities.get('description', '')}")

    # MCP / Literature integration
    if report.mcp_status:
        lines.append("")
        lines.append("MCP / Literature Integration:")
        for tool, available in report.mcp_status.items():
            status_str = (
                "[green]available[/green]" if available else "[yellow]not available[/yellow]"
            )
            lines.append(f"  {tool}: {status_str}")
        if not any(report.mcp_status.values()):
            lines.append("  (Literature search requires network access to PubMed)")

    # License status
    try:
        from .license_manager import load_license

        lic = load_license()
        lines.append("")
        lines.append("License:")
        tier_color = {
            "free": "yellow",
            "academic": "green",
            "commercial": "green",
            "trial": "yellow",
        }
        color = tier_color.get(lic.tier.value, "white")
        lines.append(f"  Tier: [{color}]{lic.tier.value}[/{color}]")
        if lic.holder:
            lines.append(f"  Holder: {lic.holder}")
        if lic.is_expired():
            lines.append("  [red]EXPIRED — renew at https://celltypepilot.io/license[/red]")
        elif lic.tier.value == "free":
            lines.append("  Upgrade: https://celltypepilot.io/license")
    except Exception:
        pass

    # Warnings
    if report.warnings:
        lines.append("")
        lines.append("Warnings:")
        for w in report.warnings:
            lines.append(f"  [WARN] {w}")

    lines.append("")
    lines.append("=" * 60)
    return "\n".join(lines)


def print_doctor():
    """Run and print the doctor report."""
    report = run_doctor()
    output = format_doctor_report(report)
    print(output)
    return report
