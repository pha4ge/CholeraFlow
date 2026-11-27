"""Interactive pipeline workflow schematic (architecture vs. results availability in this run).

Presentation-only view model. Reads the existing report context (``ctx``), the
JSON report model (``model``, see report_model.py) and the raw loaded inputs
(``raw``); never creates, changes or infers scientific results.
"""


def build_workflow_view(ctx, model, raw):
    """Return a JSON-serialisable dict consumed by the page partial(s)."""
    return {}
