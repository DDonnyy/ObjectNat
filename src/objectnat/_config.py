import sys
from typing import Literal

from iduedu import config as iduedu_config
from loguru import logger


class Config:
    """
    A configuration class to manage global settings for the application, such as Overpass API URL,
    timeouts, and logging options.

    Attributes
    ----------
    overpass_url : str
        URL for accessing the Overpass API. Defaults to "http://lz4.overpass-api.de/api/interpreter".
    timeout : int or None
        Timeout in seconds for API requests. If None, no timeout is applied.
    enable_tqdm_bar : bool
        Enables or disables progress bars (via tqdm). Defaults to True.
    logger : Logger
        Logging instance to handle application logging.

    Methods
    -------
    change_logger_lvl(lvl: Literal["TRACE", "DEBUG", "INFO", "WARN", "ERROR"])
        Changes the logging level to the specified value.
    set_overpass_url(url: str)
        Sets a new Overpass API URL.
    set_timeout(timeout: int)
        Sets the timeout for API requests.
    set_enable_tqdm(enable: bool)
        Enables or disables progress bars in the application.
    """

    def __init__(
        self,
        overpass_url="http://lz4.overpass-api.de/api/interpreter",
        timeout=None,
        enable_tqdm_bar=True,
    ):
        self.overpass_url = overpass_url
        self.timeout = timeout
        self.enable_tqdm_bar = enable_tqdm_bar
        self.logger = logger
        self.iduedu_config = iduedu_config

    def change_logger_lvl(self, lvl: Literal["TRACE", "DEBUG", "INFO", "WARN", "ERROR"]):
        self.logger.remove()
        self.logger.add(sys.stderr, level=lvl)
        self.iduedu_config.change_logger_lvl(lvl)

    def set_overpass_url(self, url: str):
        self.overpass_url = url
        self.iduedu_config.set_overpass_url(url)

    def set_timeout(self, timeout: int):
        self.timeout = timeout
        self.iduedu_config.set_timeout(timeout)

    def set_enable_tqdm(self, enable: bool):
        self.enable_tqdm_bar = enable
        self.iduedu_config.set_enable_tqdm(enable)


config = Config()
config.change_logger_lvl("INFO")
