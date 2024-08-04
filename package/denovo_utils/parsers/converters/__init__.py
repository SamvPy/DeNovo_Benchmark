"""Collection of search result parsers."""

from .base import DenovoEngineConverter
from .spectralis import SpectralisParser

__all__ = ["DenovoEngineConverter", "SpectralisParser"]
