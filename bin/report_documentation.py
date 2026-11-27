"""Methods, software versions, provenance, downloads, about and guide pages.

Presentation-only view model. Reads the existing report context (``ctx``), the
JSON report model (``model``, see report_model.py) and the raw loaded inputs
(``raw``); never creates, changes or infers scientific results.
"""


def build_documentation_view(ctx, model, raw):
    """Return a JSON-serialisable dict consumed by the page partial(s)."""
    return {}
