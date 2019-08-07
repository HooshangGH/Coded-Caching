from __future__ import division
import numpy as np
#import colorama
#from colorama import Fore
from pulp import *
import copy
import os
from itertools import combinations


clear = lambda: os.system('clear')

class UserGroup:
    def __init__(self, users):
        self.users = users
        self.key = ','.join(str(k) for k in users)
        self.time_intervals = []
        self.carriable_files = [[] for _ in range(max(users)+1)]
        self.val = 0


class X_U_l:
    def __init__(self, users, interval, val, key):
        self.users = users
        self.interval = interval
        self.val = val
        self.key = key
    def __lt__(self, other):
        if self.interval < other.interval:
            return True
        elif self.interval > other.interval:
            return False
        return len(self.users) >= len(other.users)


class Decentralized_CC:
    def __init__(self, K, N, F, n, r,
                 cache_contents,
                 requested_files,
                 arrival_times,
                 deadlines,
                 threshold_eta):
        self.inv_n = 1 / n
        self.K = K
        self.N = N
        self.F = F
        self.r = r
        self.cache_contents = cache_contents
        self.requested_files = requested_files
        self.arrivale_times = arrival_times
        self.deadlines = deadlines
        self.threshold_eta = threshold_eta

        self.Omega = [[] for k in range(K)]
        self.active_users = []
        self.current_all_times = []
        # self.active_users_intervals = []
        self.len_current_intervals = []

        self.first_constraint = []
        self.second_constraint_left = {}
        self.second_constraint_right = {}
        self.third_constraint = {}
        self.num_var_x = 0
        self.num_var_y = 0
        self.label_x = []
        self.x_U = []

        self.submitted_user_groups = {}
        self.user_benefits = [0 for _ in range(K)]


    def is_one_but_all(self, user_group):
        key_ug = ','.join(str(k) for k in user_group)
        # if key_ug in self.submitted_user_groups:
        #     return True, self.submitted_user_groups[key_ug]
        tmp_UG = UserGroup(user_group)
        for k in user_group:
            carriable_files = self.Omega[k]
            for kk in user_group:
                if kk != k:
                    carriable_files = [f for f in carriable_files if f in self.cache_contents[kk]]
            if len(carriable_files) > 0:
                tmp_UG.carriable_files[k] = carriable_files
            else:
                return False, None
        return True, tmp_UG

    def make_time_intervals(self):
        last_user = max(self.active_users)
        curr_arrival_times = [self.arrivale_times[last_user]]
        curr_deadlines = [self.deadlines[k] for k in self.active_users if self.deadlines[k] > self.arrivale_times[last_user]]
        self.current_all_times = curr_arrival_times + curr_deadlines
        self.current_all_times = list(set(self.current_all_times))
        self.current_all_times.sort()
        # self.active_users_intervals = [[-1,-1] for _ in range(last_user+1)]
        # for k in self.active_users:
        #     self.active_users_intervals[k][0] = self.current_all_times.index(self.arrivale_times[k])
        #     self.active_users_intervals[k][1] = self.current_all_times.index(self.deadlines[k])
        self.len_current_intervals = np.diff(self.current_all_times)

    def find_intervals(self, user_group):
        active_intervals = []
        for nl in range(len(self.len_current_intervals)):
            nS = self.current_all_times[nl]
            nE = self.current_all_times[nl+1]
            active_in_this_interval = [k for k in user_group if self.arrivale_times[k] <= nS and self.deadlines[k] >= nE]
            if len(active_in_this_interval) == len(user_group):
                active_intervals.append(nl)
        return active_intervals



    def make_constraints(self):
        self.num_var_x = 0
        self.num_var_y = 0
        self.first_constraint = [[] for _ in range(len(self.len_current_intervals))]
        self.second_constraint_left = {}
        self.second_constraint_right = {}
        self.third_constraint = {}
        self.all_user_groups = {}
        self.x_U = []
        self.label_x = []
        self.label_y = []
        for len_ug in range(len(self.active_users)):
            for usr_grp in combination(list(self.active_users), len_ug + 1):
                active_intervals_ug = self.find_intervals(usr_grp)
                if len(active_intervals_ug) > 0:
                    is_valid, user_group = self.is_one_but_all(usr_grp)
                    if is_valid:
                        user_group.time_intervals = active_intervals_ug
                        self.all_user_groups[user_group.key] = user_group
                        tmp_list_right = []
                        if user_group.key in self.second_constraint_right:
                            tmp_list_right = self.second_constraint_right[user_group.key]
                        for nl in range(len(active_intervals_ug)):
                            self.first_constraint[nl].append(self.num_var_x)
                            tmp_label_x = '{' + user_group.key + ' - ' + str(nl) + '}'
                            self.label_x.append(tmp_label_x)
                            tmp_list_right.append(self.num_var_x)
                            self.x_U.append(X_U_l(usr_grp, nl, 0, user_group.key))
                            self.num_var_x += 1
                        self.second_constraint_right[user_group.key] = tmp_list_right
                        for k in usr_grp:
                            sec_key = user_group.key + '-' + str(k)
                            if sec_key not in self.second_constraint_left:
                                tmp_list = []
                            else:
                                tmp_list = self.second_constraint_left[sec_key]
                            for f in user_group.carriable_files[k]:
                                tmp_list.append(self.num_var_y)
                                key_f = str(k) + '-' + str(f)
                                if key_f not in self.third_constraint:
                                    tmp_list_f = []
                                else:
                                    tmp_list_f = self.third_constraint[key_f]
                                tmp_list_f.append(self.num_var_y)
                                self.third_constraint[key_f] = tmp_list_f
                                self.label_y.append(str(k)+'-'+str(f) + '-{' + user_group.key + '}')
                                self.num_var_y += 1
                            self.second_constraint_left[sec_key] = tmp_list


    def make_n_solve_LP(self):

        self.active_users.sort()
        self.make_time_intervals()
        self.make_constraints()

        LP_Model = LpProblem('LP_DCAsynchronous', LpMinimize)
        #   LP variables
        x = LpVariable.matrix('x_var', list(range(self.num_var_x)), 0, None, LpContinuous)
        y = LpVariable.matrix('y_var', list(range(self.num_var_y)), 0, None, LpContinuous)

        # First Constraint
        for nl in range(len(self.len_current_intervals)):
            tmp_const = self.first_constraint[nl]
            if len(tmp_const) > 0:
                LP_Model += lpSum([x[i] for i in tmp_const]) <= self.len_current_intervals[nl]

        # Second Constraint:
        for key_ug in self.second_constraint_right:
            tmp_const_right = self.second_constraint_right[key_ug]
            tmp_prev_value = 0
            if key_ug in self.submitted_user_groups:
                tmp_prev_value = self.submitted_user_groups[key_ug].val
            users = self.all_user_groups[key_ug].users
            for k in users:
                sec_key = key_ug + '-' + str(k)
                tmp_const_left = self.second_constraint_left[sec_key]
                LP_Model += lpSum([y[i] for i in tmp_const_left]) <= lpSum([x[j] for j in tmp_const_right]) + tmp_prev_value
        # Third Constraint
        for key_f in self.third_constraint:
            tmp_const_third = self.third_constraint[key_f]
            LP_Model += lpSum(y[i] for i in tmp_const_third) == self.r

        # cost function
        LP_Model += lpSum(x), 'num_submitted_equations'

        # solve LP
        LP_Model.solve()
        if LpStatus[LP_Model.status] == 'Infeasible':
            return False, None

        LP_solution = []
        # print("---------")
        for i in range(self.num_var_x):
            # print("label: %s\tvlaue: %2.2f" % (self.label_x[i], value(x[i])))
            if value(x[i]) >= self.inv_n:
                self.x_U[i].val = value(x[i])
                LP_solution.append(self.x_U[i])
        # print('...\n x-solution\n')
        # for i in range(self.num_var_x):
        #     # print("label: %s\tvlaue: %2.2f" % (self.label_x[i], value(x[i])))
        #     if value(x[i]) > 0:
        #         print("label: %s\tvlaue: %2.2f" % (self.label_x[i], value(x[i])))
        # print('y-solution')
        # for j in range(self.num_var_y):
        #     if value(y[j]) > 0:
        #         print("label: %s\tvalue: %2.2f" % (self.label_y[j], value(y[j])))
        LP_solution.sort()
        return True, LP_solution

    def compute_eta(self, user_group, t):
        eta = 0
        for k in user_group:
            if self.deadlines[k] == t:
                return float('inf')
            eta += (len(self.Omega[k]) - self.user_benefits[k]) / (self.deadlines[k] - t)
        return eta

    def pick_front_load(self, LP_solution, t):
        for ns in range(len(LP_solution)):
            # print(LP_solution[ns].users)
            # print(t)
            # print(self.compute_eta(LP_solution[ns].users, t))
            # print('======================')
            if self.compute_eta(LP_solution[ns].users, t) >= self.threshold_eta and LP_solution[ns].val >= self.inv_n:
                return ns
        return -1

    def pick_front_load_modified(self, LP_solution, t):
        for ns in range(len(LP_solution)):
            # print(LP_solution[ns].users)
            # print(t)
            # print(self.compute_eta(LP_solution[ns].users, t))
            # print('======================')
            is_feasible, eta_U, modified_benefits = self.modified_eta_compute(LP_solution[ns].users,self.inv_n, t)
            if is_feasible and eta_U >= self.threshold_eta and LP_solution[ns].val >= self.inv_n:
                return ns, modified_benefits
        return -1, None

    def modified_eta_compute(self, user_group_in, incremented_value, cur_time):

        tmp_submitted_user_group = {}
        for key_ug in self.submitted_user_groups:
            tmp_submitted_user_group[key_ug] = copy.deepcopy(self.submitted_user_groups[key_ug])

        key_user_group_in = ','.join(str(k) for k in user_group_in)
        if key_user_group_in not in tmp_submitted_user_group:
            tmp_user_group = UserGroup(user_group_in)
            _, tmp_UG = self.is_one_but_all(user_group_in)
            tmp_user_group.val = incremented_value
            tmp_user_group.carriable_files = tmp_UG.carriable_files
            tmp_submitted_user_group[key_user_group_in] = tmp_user_group
        else:
            tmp_submitted_user_group[key_user_group_in].val += incremented_value

        num_var = 0
        second_const = {}
        second_val = {}
        third_const = {}
        user_vals = [[] for _ in range(self.K)]
        label_y = []
        for key_user_group in tmp_submitted_user_group:
            U_overlap = [k for k in tmp_submitted_user_group[key_user_group].users if k in user_group_in]
            for k in U_overlap:
                key_second = str(k) + '-{' + key_user_group + '}'
                tmp_second = []
                for f in tmp_submitted_user_group[key_user_group].carriable_files[k]:
                    third_key = str(k) + '-' + str(f)
                    tmp_second.append(num_var)
                    if third_key not in third_const:
                        third_const[third_key] = [num_var]
                    else:
                        third_const[third_key].append(num_var)
                    user_vals[k].append(num_var)
                    label_y.append(third_key + '-{'+key_user_group + '}')
                    num_var += 1
                second_const[key_second] = tmp_second
                second_val[key_second] = tmp_submitted_user_group[key_user_group].val

        LP_Model_eta = LpProblem('LP_Eta', LpMaximize)
        y = LpVariable.matrix('y_var', list(range(num_var)), 0, None, LpContinuous)

        # second condition
        for sec_key in second_const:
            LP_Model_eta += lpSum([y[i] for i in second_const[sec_key]]) <= second_val[sec_key]

        # third condition
        for third_key in third_const:
            LP_Model_eta += lpSum([y[i] for i in third_const[third_key]]) <= self.r

        LP_Model_eta += lpSum(y), 'benefits'

        LP_Model_eta.solve()
        # print('=======================')
        # for i in range(num_var):
        #     if value(y[i]) > 0:
        #         print('y_%s = %1.1f' % (label_y[i], value(y[i])))


        if LpStatus[LP_Model_eta.status] == 'Infeasible':
            return False, None


        eta_U = 0
        updated_benefits = [self.user_benefits[k] for k in range(self.K)]
        for k in user_group_in:
            if self.deadlines[k] <= cur_time:
                tmp_eta = 1
            else:
                tmp_eta = (len(self.Omega[k]) - self.user_benefits[k]) * self.r / (self.deadlines[k] - cur_time)
            updated_user_benefit = sum([value(y[i]) for i in user_vals[k]])
            eta_U += tmp_eta * (updated_user_benefit - self.user_benefits[k])
            updated_benefits[k] = updated_user_benefit

        return True, eta_U, updated_benefits


    def check_feasibility(self, leaving_users):
        if len(leaving_users) < 1:
            return True
        # constructing the constraints
        sec_con_left = {}
        sec_con_right = {}
        third_con = {}
        label_y = []
        num_var_y = 0
        for key_ug in self.submitted_user_groups:
            user_group = self.submitted_user_groups[key_ug].users
            overlap_users = [k for k in leaving_users if k in set(user_group)]
            if len(overlap_users) > 0:
                for k in overlap_users:
                    key_sec = str(k) + '-{' + key_ug + '}'
                    tmp_vec = []
                    for f in self.submitted_user_groups[key_ug].carriable_files[k]:
                        tmp_vec.append(num_var_y)
                        key_third = str(k) + '-' + str(f)
                        if key_third not in third_con:
                            third_con[key_third] = [num_var_y]
                        else:
                            third_con[key_third].append(num_var_y)
                        label_y.append(key_third + '{' + key_ug + '}')
                        num_var_y += 1
                    sec_con_left[key_sec] = tmp_vec
                    sec_con_right[key_sec] = self.submitted_user_groups[key_ug].val


        #   checking feasibility
        LP_Model_Feasibility = LpProblem('LP_Feasibility', LpMinimize)
        y = LpVariable.matrix('y_var', list(range(self.num_var_y)), 0, None, LpContinuous)

        # second constraint
        for key_sec in sec_con_left:
            LP_Model_Feasibility += lpSum([y[i] for i in sec_con_left[key_sec]]) <= sec_con_right[key_sec]

        # third constraint
        for key_third in third_con:
            LP_Model_Feasibility += lpSum([y[i] for i in third_con[key_third]]) == self.r

        LP_Model_Feasibility += 1, 'feasibility_check'

        #LP_Model_Feasibility.solve()
        #if LpStatus[LP_Model_Feasibility.status] == 'Infeasible':
        #    print("INFEASIBLE WHILE ALGO RETRUNS FEASIBILITY")
        #    #print()
        #    return False
        return True





    def run_algo(self):

        t_min = min(self.arrivale_times)
        t_max = max(self.deadlines) + self.inv_n
        rate = 0
        LP_sol = []

        for t in np.arange(t_min, t_max, self.inv_n):
            arrived_users = [k for k in range(self.K) if self.arrivale_times[k] == t]
            if len(arrived_users) > 0:
                for k in arrived_users:
                    d_k = (self.requested_files[k]-1)*self.F
                    self.Omega[k] = [d_k + f for f in range(self.F) if (d_k + f) not in self.cache_contents[k]]
                    self.active_users.append(k)
                is_feasible, LP_sol = self.make_n_solve_LP()
                if not is_feasible:
                    return -1, None
            idx_transmitting_ug = self.pick_front_load(LP_sol, t)
            if idx_transmitting_ug > -1:
                # print(t)
                # print(LP_sol[idx_transmitting_ug].users)
                rate += self.inv_n
                tmp_sol = LP_sol[idx_transmitting_ug]
                leaving_users = []

                if tmp_sol.key not in self.submitted_user_groups:
                    tmp_user_group = UserGroup(tmp_sol.users)
                    _, tmp_UG = self.is_one_but_all(tmp_sol.users)
                    tmp_user_group.val = 0
                    tmp_user_group.carriable_files = tmp_UG.carriable_files
                    self.submitted_user_groups[tmp_sol.key] = tmp_user_group
                self.submitted_user_groups[tmp_sol.key].val += self.inv_n
                LP_sol[idx_transmitting_ug].val -= self.inv_n
                for k in tmp_sol.users:
                    self.user_benefits[k] += self.inv_n
                    if self.user_benefits[k] >= len(self.Omega[k]) - 1e-2 and k in self.active_users:
                        self.active_users.remove(k)
                        leaving_users.append(k)
                is_leaving_users_satisfied = self.check_feasibility(leaving_users)
                if(not is_leaving_users_satisfied):
                    print("stop here!")
                    print(self.arrivale_times)
                    print(self.deadlines)
            exiting_users = [k for k in range(self.K) if self.deadlines[k] == t]
            for k in exiting_users:
                if self.user_benefits[k] < len(self.Omega[k]) - 1e-2:
                    return -1, None

        return rate, self.submitted_user_groups


    def run_algo_mod_eta(self):

        t_min = min(self.arrivale_times)
        t_max = max(self.deadlines) + self.inv_n
        rate = 0
        LP_sol = []

        for t in np.arange(t_min, t_max, self.inv_n):
            arrived_users = [k for k in range(self.K) if self.arrivale_times[k] == t]
            if len(arrived_users) > 0:
                for k in arrived_users:
                    d_k = (self.requested_files[k]-1)*self.F
                    self.Omega[k] = [d_k + f for f in range(self.F) if (d_k + f) not in self.cache_contents[k]]
                    self.active_users.append(k)
                is_feasible, LP_sol = self.make_n_solve_LP()
                if not is_feasible:
                    return -1, None
            idx_transmitting_ug, modified_benefits = self.pick_front_load_modified(LP_sol, t)
            if idx_transmitting_ug > -1:
                # print(t)
                # print(LP_sol[idx_transmitting_ug].users)
                rate += self.inv_n
                tmp_sol = LP_sol[idx_transmitting_ug]
                leaving_users = []

                if tmp_sol.key not in self.submitted_user_groups:
                    tmp_user_group = UserGroup(tmp_sol.users)
                    _, tmp_UG = self.is_one_but_all(tmp_sol.users)
                    tmp_user_group.val = 0
                    tmp_user_group.carriable_files = tmp_UG.carriable_files
                    self.submitted_user_groups[tmp_sol.key] = tmp_user_group
                self.submitted_user_groups[tmp_sol.key].val += self.inv_n
                LP_sol[idx_transmitting_ug].val -= self.inv_n
                for k in tmp_sol.users:
                    self.user_benefits[k] = modified_benefits[k]
                    if self.user_benefits[k] >= self.r * len(self.Omega[k]) - 1e-2 and k in self.active_users:
                        self.active_users.remove(k)
                        leaving_users.append(k)
                # is_leaving_users_satisfied = self.check_feasibility(leaving_users)
            exiting_users = [k for k in range(self.K) if self.deadlines[k] == t]
            for k in exiting_users:
                if self.user_benefits[k] < self.r * len(self.Omega[k]) - 1e-2:
                    return -1, None

        # Is_feasible = self.check_feasibility([_ for _ in range(self.K)])

        return rate, self.submitted_user_groups
