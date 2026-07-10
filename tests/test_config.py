"""Unit tests for the global configuration object."""

from objectnat import config
from objectnat._config import Config


def test_set_enable_tqdm_toggles_flag():
    cfg = Config()
    assert cfg.enable_tqdm_bar is True
    cfg.set_enable_tqdm(False)
    assert cfg.enable_tqdm_bar is False
    cfg.set_enable_tqdm(True)
    assert cfg.enable_tqdm_bar is True


def test_change_logger_lvl_reconfigures_without_error():
    # Swapping levels rewires the loguru sink; restore INFO afterwards so other
    # tests keep the default configured in conftest/_config.
    try:
        config.change_logger_lvl("DEBUG")
        config.change_logger_lvl("ERROR")
    finally:
        config.change_logger_lvl("INFO")
