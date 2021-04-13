from abc import ABC, abstractmethod

from utils.transport import Output


class Formulation(ABC):
    @abstractmethod
    def solve(self) -> {Output, None}:
        pass
