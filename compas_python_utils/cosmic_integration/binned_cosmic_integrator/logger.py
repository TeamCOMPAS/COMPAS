import sys
import warnings

from loguru import logger

logger.add(
    sys.stderr,
    format="|<blue>COMPAS py</>|{time:DD/MM HH:mm:ss}|{level}| {message} ",
    colorize=True,
    level="INFO",
)


warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)
