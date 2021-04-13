from formulations.formulation_abstract import Formulation
from formulations.fallback_matching import Matching
from formulations.fallback_unbounded import Unbounded
from utils.transport import Input


class Fallback(Formulation):
    fallbacks = [Unbounded, Matching]

    def __init__(self, input_data: Input):
        self.problem = None
        self.can_run = False
        for formulation in self.fallbacks:
            problem = formulation(input_data)
            if problem.can_run:
                print("Choosing fallback:", formulation.__name__)
                self.can_run = True
                self.problem = problem
                return

    def solve(self):
        if not self.can_run:
            message = "No fallback formulation applies to this instance"
            print("Warning:", message)
            return
        return self.problem.solve(), self.__class__
