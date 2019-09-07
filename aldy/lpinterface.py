# 786

# Aldy source: lpinterface.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from typing import Any, Optional, Dict, Tuple, List, Iterable, Callable, Set

import importlib
import collections

from .common import log, sorted_tuple, SOLUTION_PRECISION


def escape_name(s: str) -> str:
   return s.replace('.', '').replace('-', 'm').replace('/', '__')[:200]


class NoSolutionsError(Exception):
   pass


class Gurobi:
   """
   An abstraction aroung Gurobi Python interface.
   """

   def __init__(self, name, prev_model = None):
      self.gurobipy = importlib.import_module('gurobipy')

      self.INF = self.gurobipy.GRB.INFINITY
      self.GUROBI_STATUS = {
         getattr(self.gurobipy.GRB.status, v): v
         for v in dir(self.gurobipy.GRB.status)
         if v[:2] != '__'
      }
      if prev_model:
         self.model = prev_model
      else:
         self.model = self.gurobipy.Model(name)
         self.model.reset()


   def addConstr(self, *args, **kwargs):
      """
      Add a constraint to the model.
      """
      if 'name' in kwargs:
         kwargs['name'] = escape_name(kwargs['name'])
      c = self.model.addConstr(*args, **kwargs)
      return c


   def addVar(self, *args, **kwargs):
      """
      Add a variable to the model.

      ``vtype`` named argument stands for the variable type:

      - ``B`` for binary variable
      - ``I`` for integer variable
      - ``C`` or nothing for continuous variable.
      """
      if 'vtype' in kwargs and kwargs['vtype'] == 'B':
         kwargs['vtype'] = self.gurobipy.GRB.BINARY
      elif 'vtype' in kwargs and kwargs['vtype'] == 'I':
         kwargs['vtype'] = self.gurobipy.GRB.INTEGER
      if 'name' in kwargs:
         kwargs['name'] = escape_name(kwargs['name'])
      update = True
      if 'update' in kwargs:
         update = kwargs['update']
         del kwargs['update']
      v = self.model.addVar(*args, **kwargs)
      if update:
         self.update()
      return v


   def setObjective(self, objective, method: str = 'min'):
      """
      Sets the model objective via ``objective``.
      """
      self.objective = objective
      self.model.setObjective(
         self.objective,
         self.gurobipy.GRB.MINIMIZE if method == 'min' else self.gurobipy.GRB.MAXIMIZE)
      self.update()


   def quicksum(self, expr: Iterable):
      """
      Perform a quick summation of the iterable ``expr``.
      Much faster than Python's ``sum`` on large iterables.
      """
      return self.gurobipy.quicksum(expr)


   def update(self) -> None:
      """
      Update the model.
      Avoid calling it too much as it slows down the model construction.
      """
      self.model.update()


   def varName(self, var):
      """
      Getter for variable name.
      """
      return var.varName


   def abssum(self, vars: Iterable, coeffs: Optional[Dict[str, float]] = None):
      """
      Returns absolute sum of the ``vars``: e.g.
         :math:`\sum_i |c_i x_i|` for the set :math:`{x_1,...}`.
      where :math:`c_i` is defined in the ``coeffs`` dictionary.
      Key of the ``coeffs`` dictionary stands for the name of the variable
      accessible via ``varName`` call (1 if not defined).
      """
      vv = []
      for i, v in enumerate(vars):
         name = self.varName(v)
         coeff = 1 if coeffs is None or name not in coeffs else coeffs[name]
         absvar = self.addVar(lb=0, update=False, name=f'ABS_{name}')
         vv.append(coeff * absvar)
         self.addConstr(absvar + v >= 0, name=f'CABSL_{i}')
         self.addConstr(absvar - v >= 0, name=f'CABSR_{i}')
      self.update()
      return self.quicksum(vv)


   def prod(self, res, terms):
      """
      Sets :math:`res = \prod terms` where ``terms`` variables are all binary
      by adding appropriate linear constraints.
      Returns ``res``.
      """
      for v in terms:
         self.addConstr(res <= v)
      self.addConstr(res >= self.quicksum(terms) - (len(terms) - 1))
      return res


   def solve(self, init: Optional[Callable] = None) -> Tuple[str, float, dict]:
      """
      Solve the model. Assumes that objective is set.

      Additional parameters of the solver can be set via ``init`` function that takes
      the model instance as a sole argument.

      Returns:
         tuple[str, float]: Tuple describing the status of the solution and the objective value.

      Raises:
         :obj:`NoSolutionsError` if the model is infeasible.
      """

      self.model.params.outputFlag = 0
      self.model.params.logFile = ''
      if init is not None:
         init(self.model)
      self.model.optimize()

      status = self.GUROBI_STATUS[self.model.status]
      if self.model.status == self.gurobipy.GRB.INFEASIBLE:
         raise NoSolutionsError(status)
      return status.lower(), self.model.objVal


   def solveAll(self,
                keys: dict,
                init: Optional[Callable] = None) -> Tuple[str, float, List[tuple]]:
      """
      Solve the model. Assumes that objective is set.
      Returns the list of all combinations of the variables ``keys`` that minimize the objective.

      Additional parameters of the solver can be set via ``init`` function that takes
      the model instance as a sole argument.

      Returns:
         tuple[str, float, list[tuple[any]]]: Tuple describing the status of the solution and the objective value.
      """
      status, opt_value = self.solve(init)
      sol = sorted_tuple(set(a for a, v in keys.items() if self.getValue(v)))
      yield status, opt_value, sol
      yield from get_all_solutions(self, keys, opt_value, sol)


   def solutions(self, 
                 gap: float = 0, 
                 best_obj: Optional[float] = None, 
                 limit = None,
                 init: Optional[Callable] = None):
      def model_init(m):
         m.params.poolSearchMode = 2
         m.params.poolSolutions = limit or 2000000000
         m.params.poolGap = gap
      self.solve(init=model_init)
      for soli in range(self.model.solCount):
         self.model.params.solutionNumber = soli
         vv = {v.VarName: v
               for v in self.model.getVars()
               if v.vtype == self.gurobipy.GRB.BINARY and round(v.xn) > 0}
         yield (self.GUROBI_STATUS[self.model.status], 
                self.model.poolObjVal,
                sorted_tuple(set(vv.keys())))


   def getValue(self, var):
      """
      Get the value of the solved variable.
      Automatically adjusts the return type depending on the variable type.
      """
      if var.vtype == self.gurobipy.GRB.BINARY:
         return round(var.x) > 0
      if var.vtype == self.gurobipy.GRB.INTEGER:
         return int(round(var.x))
      else:
         return var.x


   def dump(self, file):
      self.model.write(file)


class SCIP(Gurobi):
   """
   An abstraction around SCIP `PySCIPopt` Python interface.
   """

   def __init__(self, name):
      self.pyscipopt = importlib.import_module('pyscipopt')
      self.INF = 1e20
      self.SOLVER_PRECISON = 1e-5 # Use Gurobi's default precision
      self.model = self.pyscipopt.Model(name)


   def update(self):
      pass


   def addConstr(self, *args, **kwargs):
      if 'name' in kwargs:
         kwargs['name'] = escape_name(kwargs['name'])
      return self.model.addCons(*args, **kwargs)


   def addVar(self, *args, **kwargs):
      if 'name' in kwargs:
         kwargs['name'] = escape_name(kwargs['name'])
      if 'update' in kwargs:
         del kwargs['update']
      return self.model.addVar(*args, **kwargs)


   def setObjective(self, objective, method: str = 'min'):
      self.objective = objective
      self.model.setObjective(
         self.objective,
         'minimize' if method == 'min' else 'maximize'
      )


   def quicksum(self, expr):
      return self.pyscipopt.quicksum(expr)


   def solve(self, init: Optional[Callable] = None) -> Tuple[str, float]:
      # self.model.setRealParam('limits/time', 120)
      self.model.hideOutput()
      if init is not None:
         init(self.model)
      self.model.optimize()

      status = self.model.getStatus()
      if status == 'infeasible':
         raise NoSolutionsError(status)
      return status.lower(), self.model.getObjVal()


   def varName(self, var):
      return var.name


   def getValue(self, var):
      x = self.model.getVal(var)
      if var.vtype() == 'BINARY':
         return round(x) > 0
      if var.vtype() == 'INTEGER':
         return int(round(x))
      else:
         return x


   def dump(self, file):
      self.model.writeProblem(file)

   
   def variables(self):
      return self.model.getVars() 


   def is_binary(self, v):
      return v.vtype() == 'BINARY'


   def change_model(self):
      self.model.freeTransform()


   def solutions(self, 
                 gap: float = 0, 
                 best_obj: Optional[float] = None, 
                 limit = None,
                 iteration = 0, 
                 init: Optional[Callable] = None):
      try:
         status, obj = self.solve(init)
         best_obj = obj if best_obj is None else best_obj
         if status != 'optimal':
            return 
         ub = (1 + gap) * best_obj
         if abs(obj - ub) >= self.SOLVER_PRECISON and obj > ub:
            return
         
         vv = {self.varName(v): v for v in self.variables() if self.is_binary(v) and self.getValue(v) == 1}
         yield status, obj, sorted_tuple(set(vv.keys()))
         
         if not limit or iteration + 1 < limit:
            self.change_model()
            self.addConstr(self.quicksum(vv.values()) <= len(vv) - 1)
            yield from self.solutions(gap, best_obj, limit, iteration + 1, init)
      except NoSolutionsError:
         return


class CBC(SCIP):
   """
   An abstraction around CBC via Google's `or-tools` Python interface.
   """

   def __init__(self, name):
      self.ortools = importlib.import_module('ortools.linear_solver.pywraplp')
      self.model = self.ortools.Solver(name, self.ortools.Solver.CBC_MIXED_INTEGER_PROGRAMMING)
      self.INF = self.model.infinity()
      self.SOLVER_PRECISON = 1e-5 # Use Gurobi's default precision
      self.STATUS = collections.defaultdict(lambda: 'UNKNOWN', {
         self.ortools.Solver.OPTIMAL: 'OPTIMAL',
         self.ortools.Solver.FEASIBLE: 'FEASIBLE',
         self.ortools.Solver.INFEASIBLE: 'INFEASIBLE',
         self.ortools.Solver.UNBOUNDED: 'UNBOUNDED',
         self.ortools.Solver.ABNORMAL: 'ABNORMAL',
         self.ortools.Solver.NOT_SOLVED: 'NOT_SOLVED',
      })


   def update(self):
      pass


   def addConstr(self, *args, **kwargs):
      if 'name' in kwargs:
         kwargs['name'] = escape_name(kwargs['name'])
      return self.model.Add(*args, **kwargs)


   def addVar(self, *args, **kwargs):
      name = escape_name(kwargs.get('name', ''))
      lb = kwargs.get('lb', 0)
      ub = kwargs.get('ub', self.INF)
      if 'vtype' in kwargs and kwargs['vtype'] == 'B':
         v = self.model.BoolVar(name)
      elif 'vtype' in kwargs and kwargs['vtype'] == 'I':
         v = self.model.IntVar(lb, ub, name)
      else:
         v = self.model.NumVar(lb, ub, name)
      return v


   def setObjective(self, objective, method: str = 'min'):
      self.objective = objective
      if method == 'min':
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
      if not self.model.VerifySolution(self.SOLVER_PRECISON, True):
         raise NoSolutionsError(status)
      return self.STATUS[status].lower(), self.model.Objective().Value()


   def varName(self, var):
      return var.name()


   def getValue(self, var):
      x = var.solution_value()
      if var.integer():
         x = int(round(x))
         if abs(var.lb()) < SOLUTION_PRECISION and abs(1 - var.ub()) < SOLUTION_PRECISION:
            return x > 0
         else:
            return x
      else:
         return x


   def dump(self, file):
      log.warn('Dumping not supported with CBC solver')
      pass


   def variables(self):
      return self.model.variables() 


   def is_binary(self, v):
      return isinstance(self.getValue(v), bool)


   def change_model(self):
      pass


class MIPCL(SCIP):
   """
   An abstraction around CBC via Google's `or-tools` Python interface.
   """

   def __init__(self, name):
      self.mip = importlib.import_module('mipcl_py.mipshell.mipshell')
      self.model = self.mip.Problem(name)
      self.INF = self.mip.VAR_INF
      self.SOLVER_PRECISON = 1e-5 # Use Gurobi's default precision


   def update(self):
      pass


   def addConstr(self, *args, **kwargs):
      pass


   def addVar(self, *args, **kwargs):
      if 'vtype' in kwargs:
         kwargs['type'] = kwargs['vtype']
         del kwargs['vtype']
      if 'type' in kwargs and kwargs['type'] == 'B':
         kwargs['type'] = self.mip.BIN
      elif 'type' in kwargs and kwargs['type'] == 'I':
         kwargs['type'] = self.mip.INT
      elif 'type' in kwargs:
         kwargs['type'] = self.mip.REAL
      if 'name' in kwargs:
         kwargs['name'] = escape_name(kwargs['name'])
      v = self.mip.Var(*args, **kwargs)
      return v


   def setObjective(self, objective, method: str = 'min'):
      self.objective = objective
      if method == 'min':
         self.model.minimize(self.objective)
      else:
         self.model.maximize(self.objective)


   def quicksum(self, expr):
      return self.mip.sum_(expr)


   def solve(self, init: Optional[Callable] = None) -> Tuple[str, float]:
      if init is not None:
         init(self.model)
      self.model.optimize()

      if not model.is_solutionOptimal:
         raise NoSolutionsError(status)
      return 'optimal', self.model.getObjVal()


   def varName(self, var):
      return var.name


   def getValue(self, var):
      x = var.val
      if var.type == self.mip.INT:
         return int(round(x))
      elif var.type == self.mip.BIN:
         return int(round(x)) > 0
      else:
         return x


   def dump(self, file):
      log.warn('Dumping not supported with MIPCL solver')
      pass


   def vars(self):
      return self.model.vars


   def is_binary(self, v):
      return var.type == self.mip.BIN


   def change_model(self):
      pass



def model(name: str, solver: str):
   """
   Create the ILP solver instance for a model named ``name``.
   If ``solver`` is ``'any'``, this function will attempt to use
   Gurobi, and will fall back on SCIP if Gurobi fails.

   Raises:
      :obj:`Exception` if no solver is found.
   """

   def test_gurobi(name):
      """
      Tests if Gurobi is present. Requires Gurobi 7+.
      """
      try:
         model = Gurobi(name)
         log.trace('Using Gurobi')
      except ImportError as e:
         log.warn('Gurobi not found. Please install Gurobi and gurobipy Python package.')
         log.error('{}', e)
         model = None
      return model

   def test_scip(name):
      """
      Tests if SCIP is present. Requires `PySCIPopt`.
      """
      try:
         model = SCIP(name)
         log.trace('Using SCIP')
      except ImportError as e:
         log.warn('SCIP not found. Please install SCIP and pyscipopt Python package.')
         log.error('{}', e)
         model = None
      return model

   def test_cbc(name):
      """
      Tests if OR-Tools are present. Requires Google's `ortools`.
      """
      try:
         model = CBC(name)
         log.trace('Using CBC')
      except ImportError as e:
         log.warn('CBC (Google OR-Tools) not found. Please install ortools Python package.')
         log.error('{}', e)
         model = None
      return model


   if solver == 'any':
      model = test_gurobi(name)
      if model is None:
         model = test_scip(name)
      if model is None:
         model = test_cbc(name)
      if model is None:
         raise Exception('No IP solver found. Aldy cannot solve any problems without matching IP solver. Please try installing Gurobi or SCIP.')
      return model
   else:
      fname = 'test_' + solver
      if fname in locals():
         return locals()[fname](name)
      else:
         raise Exception('IP solver {} is not supported'.format(solver))

