# 786
# Aldy source: lpinterface.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from typing import Optional, Dict, Tuple, Iterable, Callable
import importlib
import collections

from .common import log, sorted_tuple, SOLUTION_PRECISION


SOLVER_PRECISON = 1e-5
"""Default solver precision"""


def escape_name(s: str, d: Optional[dict] = None) -> str:
    """
    Escape variable names to conform with the various solver requirements.
    """
    s = s.replace(".", "").replace("-", "m").replace("#", "__").replace(">", "")[:200]
    # Ensure that all names are unique
    if d is not None:
        d[s] += 1
        if d[s] > 1:
            return s + f"_{d[s]}"
    return s


class NoSolutionsError(Exception):
    """Raised if a model is infeasible."""

    pass


class Gurobi:  # pragma: no cover
    """Wrapper around Gurobi's Python interface (:py:mod:`gurobipy`)."""

    def __init__(self, name, prev_model=None):
        self.gurobipy = importlib.import_module("gurobipy")
        self.names = collections.defaultdict(int)

        self.env = self.gurobipy.Env(empty=True)
        self.env.setParam("OutputFlag", 0)
        self.env.start()
        self.INF = self.gurobipy.GRB.INFINITY
        self.GUROBI_STATUS = {
            getattr(self.gurobipy.GRB.status, v): v
            for v in dir(self.gurobipy.GRB.status)
            if v[:2] != "__"
        }
        if prev_model:
            self.model = prev_model
        else:
            self.model = self.gurobipy.Model(name, env=self.env)
            self.model.reset()

    def addConstr(self, *args, **kwargs):
        """Add a constraint to the model."""
        if "name" in kwargs:
            kwargs["name"] = escape_name(kwargs["name"], self.names)
        c = self.model.addConstr(*args, **kwargs)
        return c

    def addVar(self, *args, **kwargs):
        """
        Add a variable to the model.

        `vtype` is the variable type:

            - `B` for binary variable
            - `I` for integer variable
            - `C` or nothing for continuous variable.
        """
        if "vtype" in kwargs and kwargs["vtype"] == "B":
            kwargs["vtype"] = self.gurobipy.GRB.BINARY
        elif "vtype" in kwargs and kwargs["vtype"] == "I":
            kwargs["vtype"] = self.gurobipy.GRB.INTEGER
        if "name" in kwargs:
            kwargs["name"] = escape_name(kwargs["name"], self.names)
        update = True
        if "update" in kwargs:
            update = kwargs["update"]
            del kwargs["update"]
        v = self.model.addVar(*args, **kwargs)
        if update:
            self.update()
        return v

    def setObjective(self, objective, method: str = "min"):
        """Set the model objective."""
        self.objective = objective
        self.model.setObjective(
            self.objective,
            self.gurobipy.GRB.MINIMIZE
            if method == "min"
            else self.gurobipy.GRB.MAXIMIZE,
        )
        self.update()

    def quicksum(self, expr: Iterable):
        """
        Perform a quick summation of the iterable `expr`.
        Much faster than Python's `sum` on large iterables.
        """
        return self.gurobipy.quicksum(expr)

    def update(self) -> None:
        """
        Update the model.
        Avoid calling it too much as it slows down the model construction.
        """
        self.model.update()

    def varName(self, var):
        """Return a variable name."""
        return var.varName

    def abssum(self, vars: Iterable, coeffs: Optional[Dict[str, float]] = None):
        r"""
        Return the absolute sum of `vars`: e.g.
           :math:`\sum_i |c_i x_i|` for the set :math:`{x_1,...}`.
        where :math:`c_i` is defined in the `coeffs` dictionary.

        Key of the `coeffs` dictionary stands for the name of the variable
        (should be accessible via `varName` call; 1 if not defined).
        """
        vv = []
        for i, v in enumerate(vars):
            name = self.varName(v)
            coeff = 1 if coeffs is None or name not in coeffs else coeffs[name]
            absvar = self.addVar(lb=0, update=False, name=f"ABS_{name}")
            vv.append(coeff * absvar)
            self.addConstr(absvar + v >= 0, name=f"CABSL_{i}")
            self.addConstr(absvar - v >= 0, name=f"CABSR_{i}")
        self.update()
        return self.quicksum(vv)

    def prod(self, res, terms):
        r"""
        Ensure that :math:`res = \prod terms`
        (where `terms` is a sequence of binary variables)
        by adding the appropriate linear constraints.
        Returns `res`.
        """
        for v in terms:
            self.addConstr(res <= v, name="PROD")
        self.addConstr(res >= self.quicksum(terms) - (len(terms) - 1), name="PROD")
        return res

    def solve(self, init: Optional[Callable] = None) -> Tuple[str, float]:
        """
        Solve the model. Assumes that the objective is set.

        Additional parameters of the solver can be set via `init` function that takes
        the model instance as the sole argument.

        :returns: Status of the solution and the objective value.
        :raise: :py:class:`NoSolutionsError` if the model is infeasible.
        """

        self.model.params.outputFlag = 0
        self.model.params.logFile = ""
        if init is not None:
            init(self.model)
        self.model.optimize()

        status = self.GUROBI_STATUS[self.model.status]
        if self.model.status == self.gurobipy.GRB.INFEASIBLE:
            raise NoSolutionsError(status)
        return status.lower(), self.model.objVal

    def getValue(self, var):
        """
        Get the value of the solved variable.
        Automatically adjusts the return type based on the variable type.
        """
        if hasattr(var, "vtype") and var.vtype == self.gurobipy.GRB.BINARY:
            return round(var.x) > 0
        if hasattr(var, "vtype") and var.vtype == self.gurobipy.GRB.INTEGER:
            return int(round(var.x))
        elif hasattr(var, "x"):
            return var.x
        else:
            return var.getValue()

    def dump(self, file):
        """Dump the model description (in LP format) to a file."""
        self.model.write(file)

    def variables(self):
        """Return the list of model variables."""
        return self.model.getVars()

    def is_binary(self, v):
        """Check if the variable is binary."""
        return v.vtype == self.gurobipy.GRB.BINARY

    def solutions(
        self,
        gap: float = 0,
        best_obj: Optional[float] = None,
        limit=None,
        iteration=0,
        init: Optional[Callable] = None,
    ):
        """
        Solve the model and returns the list of all optimal solutions.
        Assumes that the objective is set.
        Any solution whose score is less than (1 + `gap`) times
        the optimal solution score will be included.

        A solution is defined as a dictionary of set binary variables within
        the solution that are accessed
        by their name.

        Additional parameters of the solver can be set via `init` function that takes
        the model instance as the sole argument.

        This is a generic version that supports any solver.

        :yields: Status of the solution, the objective value and the solution itself.
        """

        try:
            status, obj = self.solve(init)
            best_obj = obj if best_obj is None else best_obj
            if status != "optimal":
                return
            ub = (1 + gap) * best_obj
            if abs(obj - ub) >= SOLVER_PRECISON and obj > ub:
                return

            vv = {
                self.varName(v): v
                for v in self.variables()
                if self.is_binary(v) and self.getValue(v) == 1
            }
            yield status, obj, sorted_tuple(set(vv.keys()))

            if not limit or iteration + 1 < limit:
                self.addConstr(self.quicksum(vv.values()) <= len(vv) - 1)
                yield from self.solutions(gap, best_obj, limit, iteration + 1, init)
        except NoSolutionsError:
            return


class CBC(Gurobi):
    """
    Wrapper around CBC's Python interface (Google's ortools).
    """

    def __init__(self, name):
        self.ortools = importlib.import_module("ortools.linear_solver.pywraplp")
        self.model = self.ortools.Solver(
            name, self.ortools.Solver.CBC_MIXED_INTEGER_PROGRAMMING
        )
        self.INF = self.model.infinity()
        self.STATUS = collections.defaultdict(
            lambda: "UNKNOWN",
            {
                self.ortools.Solver.OPTIMAL: "OPTIMAL",
                self.ortools.Solver.FEASIBLE: "FEASIBLE",
                self.ortools.Solver.INFEASIBLE: "INFEASIBLE",
                self.ortools.Solver.UNBOUNDED: "UNBOUNDED",
                self.ortools.Solver.ABNORMAL: "ABNORMAL",
                self.ortools.Solver.NOT_SOLVED: "NOT_SOLVED",
            },
        )
        self.names = collections.defaultdict(int)

    def update(self):
        pass

    def addConstr(self, *args, **kwargs):
        if "name" in kwargs:
            kwargs["name"] = escape_name(kwargs["name"], self.names)
        return self.model.Add(*args, **kwargs)

    def addVar(self, *_, **kwargs):
        name = escape_name(kwargs.get("name", ""), self.names)
        lb = kwargs.get("lb", 0)
        ub = kwargs.get("ub", self.INF)
        if "vtype" in kwargs and kwargs["vtype"] == "B":
            v = self.model.BoolVar(name)
        elif "vtype" in kwargs and kwargs["vtype"] == "I":
            v = self.model.IntVar(lb, ub, name)
        else:
            v = self.model.NumVar(lb, ub, name)
        return v

    def setObjective(self, objective, method: str = "min"):
        self.objective = objective
        if method == "min":
            self.model.Minimize(self.objective)
        else:
            self.model.Maximize(self.objective)

    def quicksum(self, expr):
        return self.model.Sum(expr)

    def solve(self, init: Optional[Callable] = None) -> Tuple[str, float]:
        if init is not None:
            init(self.model)
        status = self.model.Solve()

        if status == self.ortools.Solver.INFEASIBLE:
            raise NoSolutionsError(status)
        if not self.model.VerifySolution(SOLVER_PRECISON, True):
            raise NoSolutionsError(status)
        return self.STATUS[status].lower(), self.model.Objective().Value()

    def varName(self, var):
        return var.name()

    def getValue(self, var):
        x = var.solution_value()
        if hasattr(var, "integer") and var.integer():
            x = int(round(x))
            if (
                abs(var.lb()) < SOLUTION_PRECISION
                and abs(1 - var.ub()) < SOLUTION_PRECISION
            ):
                return x > 0
            else:
                return x
        else:
            return x

    def dump(self, file):
        content = self.model.ExportModelAsLpFormat(False)
        with open(file, "w") as f:
            f.write(content)

    def variables(self):
        return self.model.variables()

    def is_binary(self, v):
        return isinstance(self.getValue(v), bool)


def model(name: str, solver: str):
    """
    Create an ILP solver instance for a model named `name`.
    If `solver` is `'any'`, this function will try to use
    Gurobi, and will fall back on CBC if Gurobi is missing.

    :raise: :py:class:`Exception` if no solver is found.
    """

    def test_gurobi(name):  # pragma: no cover
        """Test if Gurobi is present. Requires Gurobi 7+."""
        try:
            model = Gurobi(name)
            log.trace("[lp] solver= gurobi")
        except ImportError:
            model = None
        return model

    def test_cbc(name):
        """Test if OR-Tools are present. Requires Google's `ortools`."""
        try:
            model = CBC(name)
            log.trace("[lp] solver= cbc")
        except ImportError:
            model = None
        return model

    if solver == "any":
        model = test_cbc(name)
        if model is None:
            model = test_gurobi(name)
        if model is None:
            raise Exception(
                "No ILP solver found. Aldy cannot operate without an ILP solver. "
                + "Please install Gurobi or Google OR Tools."
            )
        return model
    else:
        fname = "test_" + solver
        if fname in locals():
            return locals()[fname](name)
        else:
            raise Exception("ILP solver {} is not supported".format(solver))
