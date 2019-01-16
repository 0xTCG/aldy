# 786

# Aldy source: lpinterface.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from builtins import object

import importlib

from .common import log


def model(name, solver):
   def test_gurobi(name):
      try:
         model = Gurobi(name)
         log.debug('Using Gurobi')
      except ImportError as e:
         log.warn('Gurobi not found. Please install Gurobi and gurobipy Python package.')
         log.error('{}', e)
         model = None
      return model

   # CPLEX does not support our quadratic models
   # def test_cplex(name):
   #  try:
   #     model = CPLEX(name)
   #     log.warn('Using CPLEX')
   #  except:
   #     log.warn('CPLEX not found. Please install CPLEX and the following Python packages: cplex and docplex.')
   #     model = None
   #  return model

   def test_scip(name):
      try:
         model = SCIP(name)
         log.debug('Using SCIP')
      except ImportError as e:
         log.warn('SCIP not found. Please install SCIP and pyscipopt Python package.')
         log.error('{}', e)
         model = None
      return model

   if solver == 'any':
      model = test_gurobi(name)
      # if model is None:
      #  model = test_cplex(name)
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


def get_all_solutions(c, var, opt, candidates, iteration=0, mem=[]):
   if candidates in mem:
      return []
   solutions = []
   for a in candidates:
      c.changeUb(var[a], 0)
      try:
         status, obj = c.solve()
         if status == 'optimal' and abs(obj - opt) < 1e-6:
            new_solution = set(a for a, y in var.items() if round(c.getValue(y)) > 0)
            new_candidates = (set(candidates) - set([a])) | new_solution
            solutions.append(frozenset(new_solution))
            solutions += get_all_solutions(c, var, opt, new_candidates, iteration + 1, mem)
      except NoSolutionsError:
         pass
      c.changeUb(var[a], 1)
      mem.append(candidates)
   return solutions


class NoSolutionsError(Exception):
   pass


class Gurobi(object):
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
      return self.model.addConstr(*args, **kwargs)

   def addVar(self, *args, **kwargs):
      if 'vtype' in kwargs and kwargs['vtype'] == 'B':
         kwargs['vtype'] = self.gurobipy.GRB.BINARY
      v = self.model.addVar(*args, **kwargs)
      self.model.update()
      return v

   def quicksum(self, expr):
      return self.gurobipy.quicksum(expr)

   def update(self):
      self.model.update()

   def varName(self, var):
      return var.varName

   def abssum(self, vars, coeffs=None):
      total = 0
      for v in vars:
         name = self.varName(v)
         coeff = 1 if coeffs is None or name not in coeffs else coeffs[name]
         absvar = self.addVar()
         total += absvar * coeff
         self.addConstr(absvar + v >= 0)
         self.addConstr(absvar - v >= 0)
      return total

   def solve(self, objective=None, method='min'): # ret obj val
      if objective is not None:
         self.objective = objective
      self.model.setObjective(
         self.objective,
         self.gurobipy.GRB.MINIMIZE if method == 'min' else self.gurobipy.GRB.MAXIMIZE
      )
      self.model.update()

      # self.model.params.timeLimit = 60
      # self.model.params.threads = 2
      self.model.params.outputFlag = 0
      self.model.params.logFile = ''
      self.model.optimize()

      status = self.GUROBI_STATUS[self.model.status]
      if self.model.status == self.gurobipy.GRB.INFEASIBLE:
         raise NoSolutionsError(status)
      return status.lower(), self.model.objVal

   def solveAll(self, objective, var, method='min'):
      status, opt_value = self.solve(objective, method)
      solutions = [frozenset(a for a, y in var.items() if round(self.getValue(y)) > 0)]
      solutions += get_all_solutions(
         self, var, opt_value,
         candidates=solutions[0], iteration=0, mem=list()
      )
      solutions = [tuple(sorted(y for y in x)) for x in solutions]

      solutions = list(set(solutions))
      return status, opt_value, solutions

   def getValue(self, var):
      return var.x

   def changeUb(self, var, ub):
      var.ub = ub
      self.model.update()


class SCIP(Gurobi):
   def __init__(self, name):
      self.pyscipopt = importlib.import_module('pyscipopt')
      self.INF = 1e20
      self.model = self.pyscipopt.Model(name)

   def update(self):
      pass

   def addConstr(self, *args, **kwargs):
      return self.model.addCons(*args, **kwargs)

   def addVar(self, *args, **kwargs):
      return self.model.addVar(*args, **kwargs)

   def quicksum(self, expr):
      return self.pyscipopt.quicksum(expr)

   def solve(self, objective=None, method='min'):
      if objective is not None:
         self.objective = objective
      self.model.setObjective(
         self.objective,
         'minimize' if method == 'min' else 'maximize'
      )
      # self.model.params.timeLimit = 60
      self.model.setRealParam('limits/time', 120)
      self.model.hideOutput()
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


# CPLEX does not support our quadratic models.

# class CPLEX:
#  def __init__(self, name):
#     self.docplex = importlib.import_module('docplex.mp.model')
#     self.model = self.docplex.Model(name=name)
#     self.INF = self.model.infinity
#     print self.INF

#  def addConstr(self, *args, **kwargs):
#     return self.model.add_constraint(*args, **kwargs)

#  def addVar(self, *args, **kwargs):
#     vtype = 'c'
#     if 'vtype' in kwargs:
#        vtype = 'b' if kwargs['vtype'] == 'B' else 'c'
#        del kwargs['vtype']
#     if vtype == 'b':
#        v = self.model.binary_var(*args, **kwargs)
#     else:
#        v = self.model.continuous_var(*args, **kwargs)
#     return v

#  def quicksum(self, expr):
#     return self.model.sum(expr)

#  def update(self):
#     pass

#  def varName(self, var):
#     return var.varName

#  def abssum(self, vars, coeffs=None):
#     return self.quicksum([self.model.abs(v) for v in vars])
      
#  def solve(self, objective=None, method='min'): # ret obj val
#     if objective is not None:
#        self.objective = objective
#     if method == 'min':
#        self.model.minimize(self.objective)
#     else:
#        self.model.maximize(self.objective)
      
#     # self.model.params.timeLimit = 60
#     # self.model.params.threads = 2
#     # self.model.params.outputFlag = 0
#     # self.model.params.logFile = ''
      
#     self.model.export_as_lp(path='./MODEL_DBG.LP')
#     params = {}
#     try:
#        self.solution = self.model.solve(cplex_parameters=params, agent='local', log_output=False)
#     except:
#        self.model.export_as_lp(path='./MODEL_DBG.LP')
#        self.solution = None
#        exit(1)
#     if not self.solution:
#        raise NoSolutionsError(self.solution)
#     return 'optimal', self.solution.get_objective_value()

#  def solveAll(self, objective, var, method='min'):
#     status, opt_value = self.solve(objective, method)
#     solutions = [frozenset(a for a, y in var.iteritems() if round(self.getValue(y)) > 0)]
#     solutions += get_all_solutions(
#        self, var, opt_value,
#        candidates=solutions[0], iteration=0, mem=list()
#     )
#     solutions = map(lambda x: tuple(sorted(y for y in x)), solutions)

#     solutions = list(set(solutions))
#     return status, opt_value, solutions

#  def getValue(self, var):
#     return self.solution.get_value(var)

#  def changeUb(self, var, ub):
#     var.ub = ub
