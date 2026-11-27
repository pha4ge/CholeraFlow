"""Executive overview page (hero, headline metrics, key findings, analysis availability).

Presentation-only view model. Reads the existing report context (``ctx``), the
JSON report model (``model``, see report_model.py) and the raw loaded inputs
(``raw``); never creates, changes or infers scientific results.
"""


def build_overview_view(ctx, model, raw):
    """Return a JSON-serialisable dict consumed by the page partial(s)."""
    return {}
