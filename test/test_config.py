""" Testing input parser code. """

# Disable this pylint warning since it is the proper way
# to use a pytest fixture
# pylint: disable=redefined-outer-name

import os
from tempfile import NamedTemporaryFile
from unittest.mock import patch
import pytest

from nuchic.config import settings

YAML = """
run:
    events: 1000
    nucleus: 12C
    cascade: off
    folding: off
    output: lhe

parameters:
    binding_energy: 10
    fermi_momentum: 10
    config_type: MF

histograms:
    test:
        plot_range: [0, 1]
        bins: 100
    test2:
        bins: [0, 10, 50, 100, 200]
    test3:
        plot_range: [0, 1]
        bins: 100
        scale: log
"""


@pytest.fixture
def run_file():
    """ Create a global temporary run card. """
    tmp = NamedTemporaryFile(mode='w+', delete=False)
    tmp.write(YAML)
    tmp.close()
    yield tmp
    os.unlink(tmp.name)


def test_settings_init(run_file):
    """ Test settings initialization. """
    settings().load(run_file.name)
    assert not settings().cascade
    assert not settings().folding
    assert settings().nevents == 1000
    assert settings().output_format == 'lhe'


@patch('nuchic.histogram.Histogram.__init__', return_value=None)
def test_get_histograms(mock_hist, run_file):
    """ Test getting histograms. """
    settings().load(run_file.name)
    histograms = settings().get_histograms()
    assert len(histograms) == 3
    assert mock_hist.call_count == 3
