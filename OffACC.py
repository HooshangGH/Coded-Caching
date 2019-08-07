from __future__ import division
import numpy as np
from pulp import *
from itertools import combinations

class User_Group:
    def __init__(self, user_group, carrying_subfile):
        self.user_group = user_group
        self.carrying_subfiles = carrying_subfile
        self.active_time_intervals = []

'''
        OffACC:   Offline Asynchronous Coded Caching
            written by Hooshang Ghasemi, Oct. 2018
'''

class OffACC:

    def __init__(self, N, K, M, F, arrival_times, deadlines, cache_contents, Omega):
        '''
        :param N:               number of files
        :param K:               number of users
        :param M:               cache size normalized by file size
        :param F:               size of each file
        :param arrival_times:   arrival times of the users
        :param deadlines:       deadline of each user
        :param cache_contents:  cache content of each user
        '''
        self.N = N
        self.K = K
        self.M = M
        self.F = F
        self.arrival_time = arrival_times
        self.deadlines = deadlines
        self.cache_contents = cache_contents
        self.Omega = Omega
        self.starting_time_intervals = [0]*K
        self.end_time_intervals = [0]*K
        self.num_intervals = 0
        self.find_start_end_time_intervals()
        self.all_user_groups = []
        '''
            linear programming variables
        '''
        self.LP_model = LpProblem('LP_Asynchronous', LpMinimize)
        self.constraint_one = []        # first constraint sum_U x_U(l) <= |\Pi_l|
        self.constraint_two_left = []   # second constraint, left side sum y_{i,f}(U) <=
        self.constraint_two_right = []  # second constraint, right side <= sum x_U(l)
        self.constraint_three = [[] for _ in range(self.N*self.F)]      # third constraint, sum y_{i,f}(U) = r
        self.num_x_var = 0
        self.num_y_var = 0
        self.labels_x = []
        self.labels_y = []

    def find_start_end_time_intervals(self):
        all_times = self.arrival_time + self.deadlines
        all_times = list(set(all_times))
        all_times.sort()
        self.len_time_intervals = np.diff(all_times)
        self.num_intervals = len(self.len_time_intervals)
        for k in range(self.K):
            self.starting_time_intervals[k] = all_times.index(self.arrival_time[k])
            self.end_time_intervals[k] = all_times.index(self.deadlines[k])

    def find_time_intervals(self, usr_grp):
        '''
            find time intervals where all users in usr_grp are active
        :param usr_grp:
        :return: active_time_intervals
        '''
        left_time_interval = max([self.starting_time_intervals[k] for k in usr_grp])
        right_time_interval = min([self.end_time_intervals[k] for k in usr_grp])
        active_time_intervals = [l for l in range(left_time_interval, right_time_interval)]
        return active_time_intervals


    def does_make_all_but_one(self, user_group):
        '''
            check if the current user group can make an all-but-one type equation
        :param user_group: user group
        :return:
        '''
        carriable_missing_subfiles = [None]*(len(user_group))
        i = 0
        is_valid = True
        for k in user_group:
            # check if there is an missing sub-file which can be transmitted by this user group
            pos_sub_file = self.Omega[k]
            for j in user_group:
                if j != k:
                    pos_sub_file = [f for f in pos_sub_file if f in self.cache_contents[j]]
            carriable_missing_subfiles[i] = pos_sub_file
            if len(pos_sub_file) < 1:
                is_valid = False
                return carriable_missing_subfiles, is_valid
            i += 1
        return carriable_missing_subfiles, is_valid


    def make_user_groups(self):

        # make user groups of length len_ug
        all_users = [k for k in range(self.K)]
        self.constraint_one = [[] for _ in range(self.num_intervals)]
        for len_ug in range(self.K):
            for usr_grp in combination(all_users, len_ug+1):
                # check if user group usr_grp is corresponding to an all-but-one type of equation
                carryable_sub_files, is_valid = self.does_make_all_but_one(usr_grp)
                if is_valid:
                    tmp_usr_grp = User_Group(usr_grp, carryable_sub_files)
                    tmp_usr_grp.active_time_intervals = self.find_time_intervals(usr_grp)
                    if len(tmp_usr_grp.active_time_intervals) > 0:
                        label_usr_grp = ','.join(str(k) for k in usr_grp)
                        self.all_user_groups.append(tmp_usr_grp)
                        #   x_U variables corresponding to usr_grp
                        count_s = self.num_x_var
                        for l in tmp_usr_grp.active_time_intervals:
                            tmp_vec = self.constraint_one[l]
                            tmp_vec.append(self.num_x_var)
                            self.constraint_one[l] = tmp_vec
                            label_x = '{' + label_usr_grp + '}' + '-' + str(l)
                            self.labels_x.append(label_x)
                            self.num_x_var += 1
                        self.constraint_two_right.append([i_x for i_x in range(count_s, self.num_x_var)])
                        tmp_con_two = [[] for _ in range(len(usr_grp))]
                        for k in range(len(usr_grp)):
                            for f in tmp_usr_grp.carrying_subfiles[k]:
                                label_y = str(k) + str(f) + '-'
                                label_y += label_usr_grp
                                self.labels_y.append(label_y)
                                tmp_con_two[k].append(self.num_y_var)
                                self.constraint_three[f].append(self.num_y_var)
                                self.num_y_var += 1
                        self.constraint_two_left.append(tmp_con_two)
    def make_n_solve_LP(self):

        self.make_user_groups()

        # optimal rate
        opt_rate = -1

        #   LP variables
        x = LpVariable.matrix('x_var', list(range(self.num_x_var)), 0, None, LpContinuous)
        y = LpVariable.matrix('y_var', list(range(self.num_y_var)), 0, None, LpContinuous)

        # first constraint
        for l in range(self.num_intervals):
            self.LP_model += lpSum(x[i] for i in self.constraint_one[l]) <= self.len_time_intervals[l]

        # second constraint
        for nu in range(len(self.all_user_groups)):
            tmp_usr_grp = self.all_user_groups[nu]
            for k in range(len(tmp_usr_grp.user_group)):
                self.LP_model += lpSum([y[i] for i in self.constraint_two_left[nu][k]]) <= lpSum([x[j] for j in self.constraint_two_right[nu]])

        # third constraint
        for f in range(self.N*self.F):
            if len(self.constraint_three[f]) > 0:
                self.LP_model += lpSum([y[i] for i in self.constraint_three[f]]) == 1

        # cost function
        self.LP_model += lpSum(x), 'Min_Rate'

        # solve LP
        self.LP_model.solve()
        if LpStatus[self.LP_model.status] == 'Infeasible':
            return None, opt_rate

        opt_rate = 0
        opt_solution = {}
        for i in range(self.num_x_var):
            tmp_val = value(x[i])
            opt_rate += tmp_val
            if tmp_val > 0:
                opt_solution[self.labels_x[i]] = tmp_val
        return opt_solution, opt_rate
