"""Local BM25 ranking and reciprocal rank fusion for tool discovery."""

from __future__ import annotations

import math
import re
from collections import Counter, defaultdict
from dataclasses import dataclass

from dnakit.exceptions import ConfigurationError

_TOKEN = re.compile(r"[\u3400-\u4dbf\u4e00-\u9fff]+|[a-z0-9]+(?:[_-][a-z0-9]+)*", re.I)
_CJK = re.compile(r"[\u3400-\u4dbf\u4e00-\u9fff]")


@dataclass(frozen=True, slots=True)
class BM25Config:
    """Okapi BM25 parameters; positive IDF, no stemming or stop-word removal."""

    k1: float = 1.2
    b: float = 0.75

    def __post_init__(self) -> None:
        for name, value in (("k1", self.k1), ("b", self.b)):
            if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
                or (not 0 < value <= 100 if name == "k1" else not 0 <= value <= 1)
            ):
                raise ConfigurationError(
                    "BM25 k1 must be in (0, 100] and b in [0, 1], both finite numbers.",
                    code="INVALID_BM25_CONFIG",
                )


def _tokens(text: str, *, query: bool = False) -> list[str]:
    """Case-fold Latin identifiers; use CJK bigrams and document unigrams.

    Keep underscore/hyphen identifiers as well as their components. Multi-character
    CJK queries use bigrams to avoid matching solely on common single characters.
    This is the same tokenization policy as DNAKit's local knowledge retrieval.
    """
    result: list[str] = []
    for word in _TOKEN.findall(text.casefold()):
        if _CJK.match(word):
            if not query or len(word) == 1:
                result.extend(word)
            result.extend(word[index : index + 2] for index in range(len(word) - 1))
        else:
            result.append(word)
            if "_" in word or "-" in word:
                result.extend(re.split(r"[_-]", word))
    return result


def _bm25_scores(
    query: str, documents: list[tuple[str, str]], config: BM25Config
) -> dict[str, float]:
    terms = set(_tokens(query, query=True))
    if not terms or not documents:
        return {}
    frequencies = [Counter(_tokens(text)) for _, text in documents]
    document_frequency: Counter[str] = Counter()
    for counts in frequencies:
        document_frequency.update(terms.intersection(counts))
    lengths = [sum(counts.values()) for counts in frequencies]
    average = sum(lengths) / len(lengths) or 1.0
    scores: dict[str, float] = {}
    for (name, _), counts, length in zip(documents, frequencies, lengths, strict=True):
        score = math.fsum(
            math.log1p(
                (len(documents) - document_frequency[term] + 0.5) / (document_frequency[term] + 0.5)
            )
            * counts[term]
            * (config.k1 + 1)
            / (counts[term] + config.k1 * (1 - config.b + config.b * length / average))
            for term in sorted(terms.intersection(counts))
        )
        if score > 0:
            scores[name] = score
    return scores


def _rrf(rankings: tuple[list[str], ...]) -> list[str]:
    """Equal-weight RRF with k=60, using canonical names to break ties."""
    scores: defaultdict[str, float] = defaultdict(float)
    for names in rankings:
        for position, name in enumerate(names, start=1):
            scores[name] += 1.0 / (60 + position)
    return sorted(scores, key=lambda name: (-scores[name], name))
