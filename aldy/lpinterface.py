# 786

# Aldy source: lpinterface.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from typing import Any, Optional, Dict, Tuple, List, Iterable, Callable, Set

import importlib

from .common import log, sorted_tuple, SOLUTION_PRECISION


class NoSolutionsError(Exception):
   pass


class Gurobi:
   """
   An abstraction aroung Gurobi Python interface.
   """

   def __init__(self, name):
      self.gurobipy = importlib.import_module('gurobipy')

      self.INF = self.gurobipy.GRB.INFINITY
      self.GUROBI_STATUS = {
         getattr(self.gurobipy.GRB.status, v): v
         for v in dir(self.gurobipy.GRB.status)
         if v[:2] != '__'
      }
      self.model = self.gurobipy.Model(name)
      self.model.reset()


   def addConstr(self, *args, **kwargs):
      """
      Add a constraint to the model.
      """
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
         absvar = self.addVar(lb=0, update=False, name=f'ABS_{v.VarName}')
         vv.append(absvar * coeff)
         self.addConstr(absvar + v >= 0, name=f'CABSL_{i}')
         self.addConstr(absvar - v >= 0, name=f'CABSR_{i}')
      self.update()
      return self.quicksum(vv)


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


   def changeUb(self, var, ub: float) -> None:
      """
      Change the upper-bound of the variable.
      """
      var.ub = ub
      self.update()

   
   def dump(self, file):
      self.model.write(file)

class SCIP(Gurobi):
   """
   An abstraction aroung SCIP `PySCIPopt` Python interface.
   """

   def __init__(self, name):
      self.pyscipopt = importlib.import_module('pyscipopt')
      self.INF = 1e20
      self.model = self.pyscipopt.Model(name)

   def update(self):
      pass

   def addConstr(self, *args, **kwargs):
      return self.model.addCons(*args, **kwargs)

   def addVar(self, *args, **kwargs):
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
      self.model.setRealParam('limits/time', 120)
      self.model.hideOutput()
      if init is not None:
         init(self.model)
      self.model.optimize()

      status = self.model.getStatus()
      if status == 'infeasible':
         raise NoSolutionsError(status)
      return status, self.model.getObjVal()

   def varName(self, var):
      return var.name

   def getValue(self, var):
      return self.model.getVal(var)

   def changeUb(self, var, ub):
      self.model.freeTransform()
      self.model.chgVarUb(var, ub)

   def dump(self, file):
      self.model.writeProblem(file)


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

   if solver == 'any':
      model = test_gurobi(name)
      if model is None:
         model = test_scip(name)
      if model is None:
         raise Exception('No IP solver found. Aldy cannot solve any problems without matching IP solver. Please try installing Gurobi or SCIP.')
      return model
   else:
      fname = 'test_' + solver
      if fname in locals():
         return locals()[fname](name)
      else:
         raise Exception('IP solver {} is not supported'.format(solver))


def get_all_solutions(model: Gurobi, 
                      var: dict, 
                      opt: float, 
                      current_sol: tuple, 
                      iteration: int = 0, 
                      mem: Optional[Set[tuple]] = None) -> Set[tuple]:
   """
   Enumerate all possible solutions that yield ``opt`` value for the ILP defined in ``model``.

   Args:
      model (:obj:`Gurobi`): 
         ILP model instance.
      var (dict[any, :obj:`Variable`]):
         Dictionary of the variables that can form the optimal solution `current_sol`.
         Key is the variable identifier that forms the solution tuple.
      opt (float):
         The optimal solution of the model.
      current_sol (tuple[any]):
         The current solution in the pool. Used to build next solution.
         Tuple members are keys in ``var`` dictionary.
      iteration (int):
         Recursion level depth. Used to terminate the functions early.
         Max allowed depth: 10.
      mem (list[tuple[any]]):
         Memoization table to avoid recomputation of the already seen solutions.

   Note:
      Variables in ``var`` **MUST BE** binary variables in order for this to work properly!

   Returns:
      list[tuple[str]]
   """

   if mem is None:
      mem = set()
   if current_sol in mem or iteration > 10:
      return
   mem.add(current_sol)
   for a in current_sol:
      # Disable a variable
      model.changeUb(var[a], 0)
      try:
         status, obj = model.solve()
         if status == 'optimal' and abs(obj - opt) < SOLUTION_PRECISION:
            new_solution = sorted_tuple(set(vn for vn, v in var.items() if model.getValue(v)))
            yield status, obj, new_solution
            yield from get_all_solutions(model, var, opt, new_solution, iteration + 1, mem)
      except NoSolutionsError:
         pass
      # Re-enable variable
      model.changeUb(var[a], 1)
