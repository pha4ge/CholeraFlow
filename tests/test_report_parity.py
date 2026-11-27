"""Guard the comparison used for actual Nextflow report outputs."""
import pytest

from check_report_parity import check_parity


REPORT = ('<style>body{color:teal}</style><section class="cf-page" id="overview">'
          '<svg class="cf-svg"><text>2</text></svg></section>'
          '<script type="application/json" id="cf-data">{"run":"demo"}</script>'
          '<script>window.reportReady=true;</script>')


def test_design_allows_run_data_to_change_but_identical_inputs_match_figures():
    actual = REPORT.replace('"demo"', '"production"')
    assert check_parity(REPORT, actual, figures=True) == (1, 1)
    assert check_parity(REPORT, actual.replace('>2<', '>10<')) == (1, 1)


@pytest.mark.parametrize('before,after,label', [
    ('color:teal', 'color:blue', 'CSS'),
    ('reportReady=true', 'reportReady=false', 'JavaScript'),
    ('id="overview"', 'id="old-overview"', 'sections'),
    ('>2<', '>10<', 'SVG'),
])
def test_design_detects_an_old_bundle_or_changed_figure(before, after, label):
    with pytest.raises(ValueError, match=label):
        check_parity(REPORT, REPORT.replace(before, after), figures=True)


def test_missing_report_bundle_is_rejected():
    with pytest.raises(ValueError, match='Missing embedded'):
        check_parity(REPORT, '<html><p>Old report</p></html>')


def test_demo_labelling_preserves_production_design_and_figures(tmp_path):
    from test_generate_report import _full, gr, TEMPLATE
    args = _full(str(tmp_path))
    args.demo = True
    demonstration = gr.render(gr.build_context(args), TEMPLATE, str(tmp_path / 'demo.html'))
    args.demo = False
    production = gr.render(gr.build_context(args), TEMPLATE, str(tmp_path / 'production.html'))
    assert 'class="cf-sidebar-demo"' in demonstration
    assert 'class="cf-sidebar-demo"' not in production
    assert check_parity(demonstration, production, figures=True) == (25, 29)
