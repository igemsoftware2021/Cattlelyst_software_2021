"""
Author: Ruben 

Uptated by Delielena Poli 07/10/2020. Master Thesis 2020 - metabolic engineering pipeline.
Modifications to the functions RunOptknock, evaluateKO and PrintResults are indicated with
in line comments.
Updated by Sara Benito 05/04/18. Master Thesis 2018- Metabolic modelling of Actinobacteria
Script to run Optknock and Robustknock algorithms that gives the possible knock-outs which would lead 
to the maximum production of a target reaction
call_Optknock_robustknock.py calls OptAndRob function
Example:
	targets_reactions=['EX_ba']
	model=cobra.io.read_sbml_model('modelo_prueba_ba.xml') #With the added reactions
	Opt_R,Rob_R=OptAndRob(model,targets_reactions,5)
	
Script works for latest cobrapy version 14.1 after the following modifications:
	-model.solution.f replace by model.optimize().objective_value
	-The chosen solver is gurobi (gurobi7.5.2_mac64.pkg). Cplex optimizer does give a null solution at every knock out selection
	and the cplex for python is not recognized when is called 
	-tilt='yes'/'no' replace as suggested by 'True'/'False'
	-from cobra.core.Solution import Solution replaced by from cobra.core.solution import Solution
	-function solve(): 4 arguments might be included in the Solution class call instead of two: self, objective_value, status and fluxes. 
	The self.solution.x and self.solution.y assigments have been comented (those arguments are deprecated and they are not used later on)
	-In function evaluateKO, FVA variable is slightly modified (see function) and minval,maxval it is now defined as: minval=FVA.minimum[targetRxn] maxval=FVA.maximum[targetRxn] 
	
"""

import cobra
from cobra.flux_analysis.variability import flux_variability_analysis
from copy import deepcopy
from pulp import LpProblem, LpMaximize, LpMinimize, LpVariable, LpAffineExpression, LpContinuous, LpBinary,\
    LpStatusOptimal, lpSum, LpStatus
from pulp.apis import gurobi_api, cplex_api, glpk_api
import pandas as pd
import logging, sys

# before;
# from cobra.core.Solution import Solution
from cobra.core.solution import Solution  # I replaced by solution in small letters
from cobra.exceptions import Infeasible


from pipeline.scripts.analysis import get_biomass_equation

logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)

# TO DO
# Try accessing Cplex directly as alternative to OptSlope code 
# Replace 'yes'/'no' tilt by True/False       
# KO_cand is sometimes '' and sometimes None. Standardize.
# Add integer cuts
# Allow single reaction (non-list) to be added in OptAndRob


def OptAndRob(model, targets, numDel, KO_cand, verbose=False):  # add again the KO_cand=''
    ### PURPOSE
    #   Performs OptKnock and RobustKnock for a model for any number of targets
    ### EXAMPLE CALL
    #   Opt_R,Rob_R = OptAndRob(model,ExchangeReactions,1)
    ### INPUT
    #   model:    CobraPy model structure
    #   targets:  A list of CobraPy reaction IDs that are to be optimised (individually)
    #   numDel:    The number of KOs for OptKnock and RobustKnock
    ### OPTIONAL INPUT
    #   KO_cand:  A list of reaction IDs for reactions that can be knocked out. Default: all reactions in model
    #     verbose:  A boolean indicator of whether intermediate results should be printed to screen. Default: False.
    ### OUTPUT
    #   Opt_R:    OptKnock results
    #   Rob_R:    RobustKnock results
    
    # If no candidate target reactions are given set all reactions in the model as targets
    if KO_cand == '':
        KO_cand = [r.id for r in model.reactions]

    # Create output structures
    Opt_R = {}
    Rob_R = {}

    # Determine number of targets
    if type(targets) == str:
        targets = [targets]
    totalNum = len(targets)  #If single reaction added this does not work!!
    
    for i,target in enumerate(targets):
        if verbose:
            print("At number "+str(i+1)+" of "+str(totalNum)+" targets")

        #OptKnock
        okmodel = Run_OptKnock(model, target, numDel, verbose=False, integralityTol=0., KO_cand=KO_cand)
        if hasattr(okmodel.solution,'KOs'):  # I replace okmodel.solution  by optimize()
            KOs = [r.id for r in okmodel.solution.KOs]
            minval,maxval,gr = evaluateKO(model, target, KOs)  # Use FBA to double-check solution
        else:
            KOs = []
            minval = 0
            maxval = 0
            gr = 0
        Opt_R[target] = {'KOs': KOs, 'minimum': minval, 'maximum': maxval, 'growthRate': gr}

        #RobustKnock
        okmodel = Run_OptKnock(model,target,numDel,tilt=True,verbose=False,integralityTol=0.,KO_cand=KO_cand)
        if hasattr(okmodel.solution,'KOs'):
            KOs = [r.id for r in okmodel.solution.KOs]
            minval,maxval,gr = evaluateKO(model,target,KOs)  # Use FBA to double-check solution
        else:
            KOs = []
            minval = 0
            maxval = 0
            gr = 0
        Rob_R[target] = {'KOs':KOs,'minimum':minval,'maximum':maxval,'growthRate':gr}
        
    return Opt_R,Rob_R

def Run_OptKnock(model,targetRxn,numDel,objRxn='',objFraction=0.2,verbose=True,tilt=False,tiltValue=0.00001,epgap=0.0001,integralityTol=0.,KO_cand = None,numSol=1):
    """
    Runs the OptKnock algorithm on the model.

    EXAMPLE CALL
       results = Run_OptKnock(model,'EX_Succinate',3)
    ----------------------------------------------------------------------
    INPUT
        model:            CobraPy model structure
        targetRxn:        Reaction ID of the reaction that is to be optimised.
        numDel:            The maximum number of allowed KOs.
    OPTIONAL INPUT
        objRxn:         Reaction ID corresponding to the cellular objective. By default OptKnock will check the cobra model for the objective reaction.
        objFraction:    The minimal fraction of normal growth that is maintained. 
        verbose:         A boolean indicator of whether intermediate results should be printed to screen. Default: False.
        tilt:            Whether or not to tilt the objective. 'no' corresponds to the original optimistic OptKnock, 'yes' corresponds to a pessimistic version more similar to RobustKnock.
        tiltValue:        The size of the (negative!) tilt of the targetRxn to the objective function.
        epgap:            Acceptable difference between the final result and the dual result. 
        integralityTol: Maximal deviation of a integer variable from an integer.
        KO_cand:        The reaction IDs of reactions that may be knocked out. Default: all reactions
        numSol:         Number of alternative solutions that are desired. These are obtained using Integer cuts.
    OUTPUTS
       okmodel:         OptKnock model structure with solution.
       dictionary:      Nested dictionary with the solutions for number of allowed KOs
    """
    #function modified by Delielena Poli in date 20/09/2020

    #Process input
    assert targetRxn in [r.id for r in model.reactions], "targetRxn is not in the model: %r" % targetRxn
    assert type(numDel) is int, "numDel is not an integer: %r" % numDel
    assert objRxn in [r.id for r in model.reactions] or not objRxn, "targetRxn is not in the model: %r" % objRxn
    assert 0<=objFraction<=1, "The objective fraction has to be between 0 and 1. The value is: %r" % objFraction
    assert verbose==True or verbose==False, "Verbose has to be True or False." % objFraction
    assert tilt==False or tilt==True,"Tilt has to be either yes or no."
    assert epgap>=0,"epgap has to be larger or equal to 0."
    assert integralityTol>=0,"integralityTol has to be larger or equal to 0."

    #Set default objective reaction
    if not objRxn:
        objRxn = [rxn.id for rxn in model.reactions if rxn.objective_coefficient>0]
        assert len(objRxn)==1,'There should be exactly one objective reaction. The number of objective reactions is: %d' % len(objRxn) #Code does not work with multiple -> terminate.
        objRxn = objRxn[0]

    #Initialize model
    okmodel = OptKnock(model,verbose) 

    #Set minimal growth rate # modified by Delielens Poli in date 20/09/2020
    #if objFraction > 0:
    #    okmodel.model.optimize()
    #    okmodel.model.reactions.get_by_id(objRxn).lower_bound = okmodel.model.optimize().objective_value*objFraction # i changed model.solution.f by model.Solution.objective_value

    #Identify KO candidates
    if KO_cand is not None:
        KO_cand = [r for r in okmodel.model.reactions if r.id in KO_cand]

    #Prepare problem
    okmodel.prepare_optknock(targetRxn,ko_candidates=KO_cand,num_deletions=numDel,tilt=tilt,tiltValue=tiltValue) #I introduce self

    #Set parameters in solver
    # epgap
    okmodel.prob.solver.epgap = epgap
    # integrality tolerance
    okmodel.prob.solver.integralityTol = integralityTol

    #Solve problem
    okmodel.solve()

    #Print results, if verbose
    if verbose:
        # adapted to the new version of PrintResults, modified by Delielena Poli in date 20/09/2020
        dictionary = PrintResults(okmodel,targetRxn) 
    # New output returned: nested dictionary of {js: {KOs}}
    return okmodel, dictionary

#This function serves to evaluate the proposed KOs in terms of minimal and maximal production values
def evaluateKO(model,targetRxn,KO_IDs): #modified by Delielena Poli in date 20/09/2020
    """
     PURPOSE
       Predicts the minimal and maximal fluxes of the target reaction and of the biomass reaction based on KOs
     EXAMPLE CALL
       evaluateKO(model,results)
     INPUT
       model:      CobraPy model structure
         targetRxn:    Reaction ID of the reaction that is to be optimised.
         KO_IDs:        Reaction IDs of the reactions that have been knocked out.
     OUTPUT
        minval:         minimal production value
        maxval:         maximal production value
        growthrate:     growth rate

    """

    #Create temporary model - make no changes to original model
    tmpmodel = deepcopy(model)

    #Change objective to biomass again
    biomass = get_biomass_equation(tmpmodel) #added by Delielena Poli in date 20/09/2020
    logging.debug(biomass)#added by Delielena Poli in date 20/09/2020
    tmpmodel.objective = biomass#added by Delielena Poli in date 20/09/2020

    #Mutant
    tmpmodel = block_reactions(tmpmodel,KO_IDs,cumulative=True)
    fba = tmpmodel.optimize()
    growthrate = fba.fluxes[biomass.id] # I change tmpmodel.solution.f by objective_value()
    logging.debug(growthrate)

    if growthrate != None and growthrate > 9.5E-4: #KOs are not lethal
        #Minimal and maximal production during maximum growth
        #Before:
        #FVA = cb.flux_analysis.variability.flux_variability_analysis(tmpmodel,fraction_of_optimum=1,reaction_list=[tmpmodel.reactions.get_by_id(targetRxn)])
        try:
            FVA = cobra.flux_analysis.flux_variability_analysis(tmpmodel,reaction_list=[tmpmodel.reactions.get_by_id(targetRxn)],fraction_of_optimum=1) #I modify the cb and put cobra. Remove the .variability		
            #minval=FVA.targetRxn['minimum']
            #maxval=FVA.targetRxn['maximum']
            minval=FVA.minimum[targetRxn]
            maxval=FVA.maximum[targetRxn]
            #Before it was like this:
            # minval = FVA[targetRxn]['minimum']
            # maxval = FVA[targetRxn]['maximum']
        except Infeasible as i:
            print('''FVA is infeasible:''', i)
            minval = float('nan')
            maxval = float('nan')
    else: #Proposed KOs are lethal, correct solution was not found.
        minval = float('nan')
        maxval = float('nan')
        growthrate = 0
    
    return minval,maxval,growthrate

def PrintResults(okmodel,targetRxn): 
    '''
     PURPOSE
       Prints the (tilted) Optknock results as well as the results of an independent evaluation 
       of the KOs using FBA. This is required as OptKnock solutions are occasionaly incorrect due to numerical issues
     EXAMPLE CALL
       PrintResults(okmodel)
     INPUT
       okmodel:     An Optknock model object.
         targetRxn:    Reaction ID of the reaction that is to be optimised.
     OUTPUT
         Prints the (tilted) OptKnock results to screen.
         Returns 
            dict_Opt:      nested dictionary per number of allowed KOs:
                Structure:      {j:
                                    {
                                    'KO':[list of KOs], 
                                    'biomass': growth rate, 
                                    'Objective': objective rate
                                    'min_target_rate': minimal production value,
                                    'max_target_rate': maximal production value
                                    }
                                }
    '''
    # funciton modified by Delielena Poli in date 20/09/2020
    # removed print to screen of FVA results and created a
    # dictionary of results instead 

    dict_Opt = {} #added by Delielena Poli in date 20/09/2020
    if okmodel.prob.objective.value() != None:
        print('\nObjective: %6.3f' %okmodel.prob.objective.value())
        print('Biomass rate %6.3f:'%okmodel.solution.x_dict[okmodel.r_biomass])
        print("Sum of mu : %6.3f" % sum(okmodel.solution.mu))
    

        print ('\nKnockouts:flux')
        for r in okmodel.solution.KOs:
            print('%s:%6.10f' % (r.id, okmodel.solution.x_dict[r]))

        print ('\nNon-binary y''s')
        for r in okmodel.model.reactions:
            val = okmodel.solution.y_dict[r]
            if val!=0 and val!=1:
                print('%s:%6.10f' % (r.id, okmodel.solution.x_dict[r]))

        #Use FBA to double-check solution
        KOs = [r.id for r in okmodel.solution.KOs]
        minval,maxval,gr = evaluateKO(okmodel.model,targetRxn,KOs) 
        #print('\nResults of independent FBA evaluation')
        #print('Biomass rate %6.3f:'%gr)
        print('Minimal target rate %6.3f:'%minval)
        print('Maximal target flux %6.3f:'%maxval)
        
        dict_Opt['KOs'] = KOs
        #dict_Opt['biomass'] = gr
        dict_Opt['Objective'] = okmodel.prob.objective.value()
        #dict_Opt['min_target_rate'] = minval
        #dict_Opt['max_target_rate'] = maxval
    else:
        pass
    return dict_Opt

def block_reactions(model,reactions,cumulative=True):
    '''
    ### PURPOSE
    #   Blocks reactions in a cobraPy model
    ### EXAMPLE CALL
    #   block_reactions(model,reactions,cumulative=False)
    ### INPUT
    #   model:      CobraPy model structure
    #     reactions:    Reaction IDs of the reactions that are to be blocked.
    ### OPTIONAL INPUT
    #     cumulative:    Whether previously blocked reactions should remain blocked. Default: True
    ### OUTPUT
    #     model:         CobraPy model structure with reactions blocked.
    '''

    rxn_list = [r for r in model.reactions if r.id in reactions]
    
    if model._trimmed_reactions is None:
        model._trimmed_reactions = {}

    if not cumulative:
        unblock_reactions(model)

    for rxn in rxn_list:
        # Running this on an already deleted reaction will overwrite the
        # stored reaction bounds.
        if rxn.id in model._trimmed_reactions:
            continue

        old_lower_bound = rxn.lower_bound
        old_upper_bound = rxn.upper_bound
        model._trimmed_reactions[rxn.id] = (old_lower_bound,
                                            old_upper_bound)
        rxn.lower_bound = 0.
        rxn.upper_bound = 0.
        model._trimmed = True

    return model

def unblock_reactions(model):
    '''
    ### PURPOSE
    #   Unblocks reactions in a cobraPy model
    ### EXAMPLE CALL
    #   unblock_reactions(model)
    ### INPUT
    #   model:      CobraPy model structure
    ### OPTIONAL INPUT
    ### OUTPUT
    #    model:         CobraPy model structure with all reactions unblocked.
    '''

    if model._trimmed_reactions is not None:
        for reaction_ID, (lower_bound, upper_bound) in \
                model._trimmed_reactions.items():
            model.reactions.get_by_id(reaction_ID).lower_bound = lower_bound
            model.reactions.get_by_id(reaction_ID).upper_bound = upper_bound

    model._trimmed_reactions = {}
    model._trimmed = False

    return model

class OptKnock(object):
    ''''
    #Code originally from OptSlope (Elad Noor)
    #Modifications:
    # 1) Changed use_glpk to solver_choice
    # 2) Removed OptSlope
    # 3) Removed Solve_FVA
    # 4) Removed prepare_FBA
    # 5) Removed functions to print results
    # 6) Removed get_reaction_by_id function - only used once, cobra_model already has a get_by_id function.
    # 7) Added objective tilting
    '''

    def __init__(self, model, verbose=False):
        self.model = deepcopy(model)
        self.verbose = verbose

        # locate the biomass reaction
        biomass_reactions = [r for r in self.model.reactions
                             if r.objective_coefficient != 0]
        if len(biomass_reactions) != 1:
            raise Exception('There should be only one single biomass reaction')
        self.r_biomass = biomass_reactions[0]
        
        self.has_flux_as_variables = False

    def prepare_optknock(self, target_reaction_id, ko_candidates=None, 
                         num_deletions=5, solver_choice='gurobi',tilt=False,tiltValue=0.00001): # solver choice replaced by gurobi
        ''' find the target reaction'''
        self.r_target = self.model.reactions.get_by_id(target_reaction_id)

        self.create_prob(sense=LpMaximize, solver_choice=solver_choice)
        self.add_primal_variables_and_constraints()
        self.add_dual_variables_and_constraints()

        if tilt==True:
            S_times_lambda = LpAffineExpression([(self.var_lambda[m], coeff)
                                             for m, coeff in self.r_target._metabolites.items()
                                             #for m, coeff in self.r_target._metabolites.iteritems()
                                             if coeff != 0])
            row_sum = S_times_lambda + self.var_w_U[self.r_target] - self.var_w_L[self.r_target]
            self.prob.constraints['dual_%s' % self.r_target.id] = (row_sum == -tiltValue)

        self.add_optknock_variables_and_constraints(tilt=tilt,tiltValue=tiltValue)

        # add the objective of maximizing the flux in the target reaction
        self.prob.setObjective(self.var_v[self.r_target])

        self.add_knockout_bounds(ko_candidates, num_deletions)
    
    def create_prob(self, sense=LpMaximize, solver_choice='gurobi'): # I have changed the solver choice by gurobi
        ''' create the LP'''
        self.prob = LpProblem('OptKnock', sense=sense)
        if solver_choice=='glpk':
            self.prob.solver = glpk_api # There are different gplk classes, choose your version
        elif solver_choice=='gurobi':
            self.prob.solver = gurobi_api.GUROBI() #Changed the path name for PuLP 2.1
        elif solver_choice=='cplex':
            self.prob.solver = cplex_api # There are different Cplex classes, choose your version
        else:
            raise Exception('%s solver not among available options' % solver_choice)
        if not self.prob.solver.available():
            raise Exception('%s not available' % solver_choice)    

    def add_primal_variables_and_constraints(self):
        ''' create the continuous flux variables (can be positive or negative)'''
        self.var_v = {}
        for r in self.model.reactions:
            self.var_v[r] = LpVariable("v_%s" % r.id,
                                       lowBound=r.lower_bound,
                                       upBound=r.upper_bound,
                                       cat=LpContinuous)

        # this flag will be used later to know if to expect the flux
        # variables to exist
        self.has_flux_as_variables = True
        
        # add the mass-balance constraints to each of the metabolites (S*v = 0)
        for m in self.model.metabolites:
            S_times_v = LpAffineExpression([(self.var_v[r], r.get_coefficient(m))
                                            for r in m.reactions])
            self.prob.addConstraint(S_times_v == 0, 'mass_balance_%s' % m.id)
    
    def add_dual_variables_and_constraints(self,M=1000):
        ''' create dual variables associated with stoichiometric constraints'''
        self.var_lambda = dict([(m, LpVariable("lambda_%s" % m.id, 
                                               lowBound=-M,
                                               upBound=M,
                                               cat=LpContinuous))
                                for m in self.model.metabolites])

        # create dual variables associated with the constraints on the primal fluxes
        self.var_w_U = dict([(r, LpVariable("w_U_%s" % r.id, lowBound=0, upBound=M, cat=LpContinuous))
                             for r in self.model.reactions])
        self.var_w_L = dict([(r, LpVariable("w_L_%s" % r.id, lowBound=0, upBound=M, cat=LpContinuous))
                             for r in self.model.reactions])

        # add the dual constraints:
        #   S'*lambda + w_U - w_L = c_biomass
        for r in self.model.reactions:
            S_times_lambda = LpAffineExpression([(self.var_lambda[m], coeff)
                                                 for m, coeff in r._metabolites.items()
                                                 #for m, coeff in r._metabolites.iteritems()
                                                 if coeff != 0])
            row_sum = S_times_lambda + self.var_w_U[r] - self.var_w_L[r]
            self.prob.addConstraint(row_sum == r.objective_coefficient, 'dual_%s' % r.id)
    
    def add_optknock_variables_and_constraints(self,tilt=False,tiltValue=0.00001,M=1000):
        ''' create the binary variables indicating which reactions knocked out'''
        self.var_y = dict([(r, LpVariable("y_%s" % r.id, cat=LpBinary))
                           for r in self.model.reactions])

        # create dual variables associated with the constraints on the primal fluxes
        self.var_mu = dict([(r, LpVariable("mu_%s" % r.id, cat=LpContinuous))
                             for r in self.model.reactions])

        # equate the objectives of the primal and the dual of the inner problem
        # to force its optimization:
        #   sum_j mu_j - v_biomass = 0
        if tilt==True:
            constr = (lpSum(self.var_mu.values()) - self.var_v[self.r_biomass] +  tiltValue*self.var_v[self.r_target] == 0)
        else:
            constr = (lpSum(self.var_mu.values()) - self.var_v[self.r_biomass] == 0)
        self.prob.addConstraint(constr, 'dual_equals_primal')

        # add the knockout constraints (when y_j = 0, v_j has to be 0)
        for r in self.model.reactions:
            # L_jj * y_j <= v_j
            self.prob.addConstraint(r.lower_bound * self.var_y[r] <= self.var_v[r], 'v_lower_%s' % r.id)
            # v_j <= U_jj * y_j
            self.prob.addConstraint(self.var_v[r] <= r.upper_bound * self.var_y[r], 'v_upper_%s' % r.id)
            
        # set the constraints on the auxiliary variables (mu):
        #    mu_j == y_j * (U_jj * w_u_j - L_jj * w_l_j)
        for r in self.model.reactions:
            w_sum = LpAffineExpression([(self.var_w_U[r], r.upper_bound),
                                        (self.var_w_L[r], -r.lower_bound)])

            # mu_j + M*y_j >= 0
            self.prob.addConstraint(self.var_mu[r] + M*self.var_y[r] >= 0, 'aux_1_%s' % r.id)
            # -mu_j + M*y_j >= 0
            self.prob.addConstraint(-self.var_mu[r] + M*self.var_y[r] >= 0, 'aux_2_%s' % r.id)
            # mu_j - (U_jj * w_u_j - L_jj * w_l_j) + M*(1-y_j) >= 0
            self.prob.addConstraint(self.var_mu[r] - w_sum + M*(1-self.var_y[r]) >= 0, 'aux_3_%s' % r.id)
            # -mu_j + (U_jj * w_u_j - L_jj * w_l_j) + M*(1-y_j) >= 0
            self.prob.addConstraint(-self.var_mu[r] + w_sum + M*(1-self.var_y[r]) >= 0, 'aux_4_%s' % r.id)

    def add_knockout_bounds(self, ko_candidates=None, num_deletions=5):
        """ 
            construct the list of KO candidates and add a constraint that
            only K (num_deletians) of them can have a y_j = 0
        """
        ko_candidate_sum_y = []
        
        if ko_candidates is None:
            ko_candidates = [r for r in self.model.reactions if r != self.r_biomass]

        for r in set(self.model.reactions).difference(ko_candidates):
            # if 'r' is not a candidate constrain it to be 'active'
            # i.e.   y_j == 1
            self.prob.addConstraint(self.var_y[r] == 1, 'active_%s' % r.id)

        # set the upper bound on the number of knockouts (K)
        #   sum (1 - y_j) <= K

        ko_candidate_sum_y = [(self.var_y[r], 1) for r in ko_candidates]
        #constr = (LpAffineExpression(ko_candidate_sum_y) >= len(ko_candidate_sum_y) - num_deletions)
        constr = (LpAffineExpression(ko_candidate_sum_y) >= len(ko_candidate_sum_y) - num_deletions)
        self.prob.addConstraint(constr, 'number_of_deletions')

    def solve(self):
        self.prob.solve()
		
        if self.prob.status != LpStatusOptimal:
            if self.verbose:
                print("LP was not solved because: ", LpStatus[self.prob.status])
            self.solution = Solution(None, None, None, None)  #It needs 4 arguments as well
        else:
			
            self.solution = Solution(self.prob.objective.value(),self.prob.status,self.var_v)# It needs 4 arguments: self, objective_value,status, fluxes  (fluxes is a panda.Series. Constains the reaction fluxes (primal values of variables)
            if self.has_flux_as_variables:
                      
                # self.solution.x = [self.var_v[r].varValue  for r in self.model.reactions] #deprecated attributes 
                # self.solution.y  = [self.var_y[r].varValue  for r in self.model.reactions]  
                self.solution.mu = [self.var_mu[r].varValue for r in self.model.reactions]

                self.solution.x_dict = {}    # x_dict deprecated attribute --> fluxes
                self.solution.y_dict = {}    # y_dict deprecated attribute --> reduced_costs Contains reaction reduced costs (dual values of variables)
                self.solution.mu_dict = {}
                for r in self.model.reactions:
                    self.solution.x_dict[r] = self.var_v[r].varValue
                    
                    self.solution.y_dict[r] = self.var_y[r].varValue
                    self.solution.mu_dict[r] = self.var_mu[r].varValue

                #Store KOs for convenience
                self.solution.KOs = [r for r in self.model.reactions if self.solution.y_dict[r]==0]

        self.solution.status = self.prob.status
        return self.solution

