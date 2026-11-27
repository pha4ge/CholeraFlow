"""Optional offline country overview for recorded sample metadata.

Maps recorded, exactly matching country names to a locally embedded,
public-domain boundary dataset. No coordinates, movement or incidence are
inferred; counts are sampled isolates (submitted records).
"""


def build_geo_view(ctx, model, raw):
    """Return {'status': 'not_available', ...} or the map view model."""
    return {'status': 'not_available', 'reason': 'Country map not configured.'}
