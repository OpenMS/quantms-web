# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

OpenMS Streamlit Template is a web application framework for building mass spectrometry (MS) analysis workflows using OpenMS/pyOpenMS. It supports both simple pyOpenMS workflows and complex multi-tool pipelines using OpenMS TOPP (The OpenMS Proteomics Pipeline) tools.

## Common Commands

```bash
# Run the app locally
streamlit run app.py

# Run tests
python -m pytest test_gui.py tests/

# Build and run with Docker (includes OpenMS TOPP tools)
docker-compose up -d --build

# Clean up old workspaces (removes workspaces older than 7 days)
python clean-up-workspaces.py
```

Note: Local runs have limited functionality. Features requiring OpenMS TOPP tools only work out of the box with Docker or when OpenMS Command Line Tools are installed separately.

## Architecture

### Core Framework (`src/workflow/`)

The workflow system is built around `WorkflowManager` as the base class with these components:

- **WorkflowManager**: Base class that orchestrates file management, parameters, command execution, and UI. Custom workflows inherit from this and override `upload()`, `configure()`, `execution()`, and `results()` methods.
- **FileManager**: Handles input/output file organization in `workflow_dir/input-files/{key}/` and `workflow_dir/results/`
- **ParameterManager**: Manages TOPP tool parameters (XML .ini files) and JSON parameters
- **CommandExecutor**: Runs external commands (TOPP tools) with threading for parallelization
- **StreamlitUI**: Provides Streamlit widgets including `upload_widget()` and `input_TOPP()` for TOPP parameter UIs
- **Logger**: Multi-level logging (minimal, commands, all) to `workflow_dir/logs/`

### Page Structure

- **Entry point**: `app.py` defines multi-page navigation using `st.navigation()`
- **Pages**: Each file in `content/` is a page that calls `page_setup()` from `src/common/common.py`, then instantiates a workflow class
- **Utility pages** (`content/`): digest.py, fragmentation.py, isotope_pattern_generator.py, peptide_mz_calculator.py provide standalone analysis tools

### Workflow Data Flow

1. `page_setup()` initializes workspace and loads parameters from `workspace/params.json`
2. Workflow class (e.g., `WorkflowTest`) inherits from `WorkflowManager`
3. Each page calls the appropriate method: `show_file_upload_section()`, `show_parameter_section()`, `show_execution_section()`, `show_results_section()`
4. Workflow execution runs in a multiprocessing.Process to avoid blocking Streamlit UI updates

### Key Patterns

- **Workspace isolation**: Each user session gets a unique workspace directory for files and parameters
- **Streamlit fragments**: Use `@st.fragment` decorator for interactive UI updates without full page reloads
- **TOPP tool execution**: `executor.run_topp("ToolName", {inputs/outputs}, {extra_params})` handles parameter files and command construction

## Configuration Files

- `settings.json`: App name, version, analytics, workspace settings
- `default-parameters.json`: Workflow default parameters
- `.streamlit/config.toml`: Streamlit server config (port 8501, 1000MB upload limit)

## Creating New Workflows

Inherit from `WorkflowManager` and implement the four core methods:

```python
from src.workflow.WorkflowManager import WorkflowManager

class MyWorkflow(WorkflowManager):
    def __init__(self):
        super().__init__("My Workflow", st.session_state["workspace"])

    def upload(self):
        self.ui.upload_widget(key="input-files", name="Input", file_types="mzML", fallback=[...])

    def configure(self):
        self.ui.input_TOPP("ToolName", custom_defaults={...}, include_parameters=[...])

    def execution(self):
        self.executor.run_topp("ToolName", {"in": [...], "out": [...]}, {...})

    def results(self):
        st.dataframe(pd.read_csv(...))
```

## Key Dependencies

- **pyOpenMS 3.5.0+**: Python bindings for OpenMS
- **Streamlit 1.43.0**: Web UI framework
- **Plotly + streamlit_plotly_events**: Interactive visualizations
- **OpenMS TOPP tools**: External command-line tools (Docker or separate install required)
