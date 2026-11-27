"""Report colour tokens: one category-to-colour mapping per report.

Colours identify supplied categories only; they never encode a new call. The
categorical order was checked with a colour-vision-deficiency validator in both
light and dark modes (adjacent-pair CVD and normal-vision separation, contrast
against the chart surface). Categorical hues are assigned in a fixed order to the
most frequent categories; further categories fold into one neutral "Other" tone
rather than generated hues, and every figure exposes the exact value in its
legend, tooltip and table. Status colours are reserved for QC-style states and
always ship with an icon and a text label in the templates.
"""
from collections import Counter

# Fixed categorical order (light surface #ffffff / dark surface #0f1c2e).
CATEGORICAL = [
    ('c1', '#2563eb', '#4f8ff7'),  # blue
    ('c2', '#e8692a', '#e0763a'),  # orange
    ('c3', '#0d9488', '#26a69a'),  # teal
    ('c4', '#b98000', '#c08d12'),  # ochre
    ('c5', '#7c3aed', '#9a7cf2'),  # violet
    ('c6', '#e11d48', '#ec5577'),  # rose
    ('c7', '#c026d3', '#cf5ade'),  # fuchsia
    ('c8', '#4d7c0f', '#6fa324'),  # olive
]
OTHER = ('other', '#8b98aa', '#6c7a90')
NOT_AVAILABLE = ('na', '#cbd5df', '#3a4b61')
AMBIGUOUS = ('amb', '#778899', '#8d9db0')

# Ordinal ramp (collection year): violet, light-surface steps all clear 2:1.
ORDINAL = [
    ('o1', '#b89dd9', '#6f558c'), ('o2', '#a582cd', '#8865ad'),
    ('o3', '#9169bc', '#a078cb'), ('o4', '#7b54a5', '#b68ee3'),
    ('o5', '#654388', '#cba7f5'), ('o6', '#4e3468', '#dfc3ff'),
]

# Sequential ramps for distance heatmaps: lightest = nearest zero (light mode).
# In dark mode the order runs from the surface towards bright values.
SEQUENTIAL = {
    'snp': {
        'light': ['#e2f5f2', '#b9e7e1', '#8ad4cd', '#55bdb8', '#23a3a5',
                  '#00898f', '#006e7a', '#035561', '#123c45'],
        'dark': ['#153b44', '#00515d', '#006874', '#007e89', '#00959c',
                 '#2aaaac', '#5dbebb', '#8ad1ca', '#b7e2db'],
    },
    'cgmlst': {
        'light': ['#e6eefc', '#c3d6fb', '#9dbaf8', '#7b9ff2', '#6187e6',
                  '#5070c9', '#4359a6', '#38457e', '#2c3255'],
        'dark': ['#2d3353', '#374378', '#41549b', '#4b68bb', '#587dd5',
                 '#6995e9', '#80acf5', '#9cc4fc', '#bedafe'],
    },
}

# Reserved status tones (fill colours; text tones live in CSS tokens).
STATUS = {
    'pass': ('s-pass', '#059669', '#34d399'),
    'warn': ('s-warn', '#d97706', '#fbbf24'),
    'fail': ('s-fail', '#dc2626', '#f87171'),
    'not_assessed': ('s-na', '#94a3b8', '#64748b'),
    'missing': ('s-missing', '#e2e8f0', '#334155'),
    'review': ('s-review', '#4f46e5', '#a5b4fc'),
}

NOT_AVAILABLE_LABEL = 'Not available'
AMBIGUOUS_LABEL = 'Ambiguous'
OTHER_LABEL = 'Other'


def entry(token, light, dark, value, count=None):
    item = {'value': value, 'color': light, 'dark': dark, 'token': token}
    if count is not None:
        item['count'] = count
    return item


def categorical_map(values):
    """Map supplied category labels to fixed tokens by frequency, then label.

    ``values`` is an iterable of labels (one per sample). The mapping depends
    only on the cohort, never on a display filter, so colours follow the
    category rather than its rank in a filtered view.
    """
    counts = Counter(value for value in values
                     if value not in (NOT_AVAILABLE_LABEL, AMBIGUOUS_LABEL))
    ordered = sorted(counts, key=lambda value: (-counts[value], value.casefold(), value))
    mapping = {}
    for index, value in enumerate(ordered):
        if index < len(CATEGORICAL):
            token, light, dark = CATEGORICAL[index]
        else:
            token, light, dark = OTHER
        mapping[value] = entry(token, light, dark, value)
    mapping[NOT_AVAILABLE_LABEL] = entry(*NOT_AVAILABLE, NOT_AVAILABLE_LABEL)
    mapping[AMBIGUOUS_LABEL] = entry(*AMBIGUOUS, AMBIGUOUS_LABEL)
    return mapping


def ordinal_map(values):
    """Map ordered labels (for example collection years) onto the ordinal ramp."""
    ordered = sorted({value for value in values
                      if value not in (NOT_AVAILABLE_LABEL, AMBIGUOUS_LABEL)},
                     key=lambda value: (len(value), value))
    mapping = {}
    steps = len(ORDINAL)
    for index, value in enumerate(ordered):
        if len(ordered) == 1:
            position = steps // 2
        else:
            position = round(index * (steps - 1) / (len(ordered) - 1))
        token, light, dark = ORDINAL[position]
        mapping[value] = entry(token, light, dark, value)
    mapping[NOT_AVAILABLE_LABEL] = entry(*NOT_AVAILABLE, NOT_AVAILABLE_LABEL)
    mapping[AMBIGUOUS_LABEL] = entry(*AMBIGUOUS, AMBIGUOUS_LABEL)
    return mapping


def status_entry(status, label=None):
    token, light, dark = STATUS.get(status, STATUS['not_assessed'])
    return entry(token, light, dark, label or status)


def legend(mapping, counts):
    """Legend rows for values present in ``counts`` (value -> sample count)."""
    rows = []
    for value, count in sorted(counts.items(), key=lambda item: (
            item[0] in (NOT_AVAILABLE_LABEL, AMBIGUOUS_LABEL), item[0].casefold(), item[0])):
        known = mapping.get(value) or entry(*OTHER, value)
        rows.append(dict(known, value=value, count=count))
    return rows


def sequential_step(value, maximum, steps=9):
    """Index of the sequential step for a finite non-negative value."""
    if maximum <= 0:
        return 0
    return min(steps - 1, int(round((value / maximum) * (steps - 1))))


def css_tokens():
    """CSS custom properties for every token, per theme (used by the template)."""
    light, dark = [], []
    for token, light_hex, dark_hex in CATEGORICAL + [OTHER, NOT_AVAILABLE, AMBIGUOUS] + ORDINAL:
        light.append('--cf-{}: {};'.format(token, light_hex))
        dark.append('--cf-{}: {};'.format(token, dark_hex))
    for token, light_hex, dark_hex in STATUS.values():
        light.append('--cf-{}: {};'.format(token, light_hex))
        dark.append('--cf-{}: {};'.format(token, dark_hex))
    for name, ramps in SEQUENTIAL.items():
        for index, (light_hex, dark_hex) in enumerate(zip(ramps['light'], ramps['dark'])):
            light.append('--cf-{}-{}: {};'.format(name, index, light_hex))
            dark.append('--cf-{}-{}: {};'.format(name, index, dark_hex))
    return {'light': ' '.join(light), 'dark': ' '.join(dark)}


def tone_tokens():
    """All tone token names, for generating `.cf-tone-*` CSS classes."""
    tokens = [token for token, _, _ in CATEGORICAL + [OTHER, NOT_AVAILABLE, AMBIGUOUS] + ORDINAL]
    tokens += [token for token, _, _ in STATUS.values()]
    return tokens
