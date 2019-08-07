import numpy as np
from pulp import *
from itertools import combinations

class User_Group:

    def __init__(self, user_group, carrying_subfile):
        self.user_group = user_group
        self.carrying_subfiles = carrying_subfile
        self.active_time_intervals = []

class XU:
    def __init__(self, user_group, interval, val):
        self.user_group = user_group
        self.interval = interval
        self.val = val
    def __lt__(self, other):
        if self.interval < other.interval:
            return True
        elif self.interval == other.interval:
            return len(self.user_group) >= len(other.user_group)
        return False


class OnACC:

    def __init__(self, K, N, F, n, eta_threshold):
        '''

        :param K:   number of users
        :param N:   number of files
        :param F:   number of subfiles
        '''
        self.K = K
        self.N = N
        self.F = F
        self.n = n
        self.eta_threshold = eta_threshold
        self.submitted_usr_grp = {}
        self.active_users = set()
        self.arrivale_times = [-1 for _ in range(K)]
        self.deadlines = [-1 for _ in range(K)]
        self.cache_contents = [[] for _ in range(K)]
        self.requested_files = [_ for _ in range(K)]
        self.Omega = [[] for _ in range(K)]

        self.start_interval = [-1 for _ in range(K)]
        self.end_interval = [-1 for _ in range(K)]
        self.len_time_intervals = []

    def remove_user(self, usr):
        if usr in self.active_users:
            self.active_users.remove(usr)

    def update_submitted_usr_grp(self, key_usr_grp, added_val_usr_grp):
        if key_usr_grp not in self.submitted_usr_grp:
            self.submitted_usr_grp[key_usr_grp] = 0
        self.submitted_usr_grp[key_usr_grp] += added_val_usr_grp

    def make_intervals(self):
        current_arrival = [self.arrivale_times[k] for k in self.active_users]
        current_deadline = [self.deadlines[k] for k in self.active_users]
        self.all_times = current_arrival + current_deadline
        min_time = min([self.arrivale_times[k] for k in self.active_users])
        self.all_times = [t for t in self.all_times if t >= min_time]
        self.all_times = list(set(self.all_times))
        self.all_times.sort()
        for k in range(self.K):
            if k in self.active_users:
                self.start_interval[k] = self.all_times.index(self.arrivale_times[k])
                self.end_interval[k] = self.all_times.index(self.deadlines[k])
            else:
                self.start_interval[k] = -1
                self.end_interval[k] = -1
        self.len_time_intervals = np.diff(self.all_times)

    def does_make_all_but_one(self, user_group):
        carriable_sub_files = [[] for _ in range(self.K)]
        for k in user_group:
            tmp_subfiles = self.Omega[k]
            for n in user_group:
                if n != k:
                    tmp_subfiles = [f for f in tmp_subfiles if f in self.cache_contents[n]]
            if len(tmp_subfiles) < 1:
                return carriable_sub_files, False
            carriable_sub_files[k] = tmp_subfiles
        return carriable_sub_files, True

    def find_active_intervals(self, user_group):
        active_intervals = []
        for l in range(len(self.len_time_intervals)):
            nS = self.all_times[l]
            nE = self.all_times[l+1]
            nI = len([k for k in user_group if self.arrivale_times[k] <= nS and self.deadlines[k] >= nE])
            if nI == len(user_group):
                active_intervals.append(l)
        return active_intervals

    def build_user_groups(self):

        self.constraint_one = [[] for _ in range(len(self.len_time_intervals))]
        self.constraint_two_left = []
        self.constraint_two_right = []
        self.constraint_three = [[] for _ in range(self.N * self.F)]
        self.label_x = []
        self.label_usr_grp = []
        self.U_x = []
        self.l_x = []
        self.num_var_x = 0
        self.num_var_y = 0

        self.current_interval = self.all_times.index(max([self.arrivale_times[k] for k in self.active_users]))
        for len_usr_grp in range(len(self.active_users)):
            for usr_grp in combination(list(self.active_users), len_usr_grp+1):
                active_intervals = self.find_active_intervals(usr_grp)
                tmp_x_U_l = []
                tmp_y_i_f_U = [[] for _ in range(self.K)]
                if len(active_intervals) > 0:
                    carriable_subfiles, does_it = self.does_make_all_but_one(usr_grp)
                    if does_it:
                        label_usr_grp = ','.join(str(k) for k in usr_grp)
                        tmp_usr_grp = User_Group(usr_grp, carriable_subfiles)
                        tmp_usr_grp.active_time_intervals = active_intervals
                        for nl in active_intervals:
                            if nl >= self.current_interval:
                                self.constraint_one[nl].append(self.num_var_x)
                                label_x = '{' + label_usr_grp + '}' + '-' + str(nl)
                                tmp_x_U_l.append(self.num_var_x)
                                self.label_x.append(label_x)
                                self.U_x.append(usr_grp)
                                self.l_x.append(nl)
                                self.num_var_x += 1
                        self.constraint_two_right.append(tmp_x_U_l)
                        for k in usr_grp:
                            for f in carriable_subfiles[k]:
                                tmp_y_i_f_U[k].append(self.num_var_y)
                                self.constraint_three[f].append(self.num_var_y)
                                self.num_var_y += 1
                        self.constraint_two_left.append(tmp_y_i_f_U)
                        self.label_usr_grp.append(label_usr_grp)


    def user_arrived(self, usr, arrival_time, deadline, req_file, cache_content):
        self.active_users.add(usr)
        self.arrivale_times[usr] = arrival_time
        self.deadlines[usr] = deadline
        self.requested_files[usr] = req_file
        self.cache_contents[usr] = cache_content
        req_sub_files = [((req_file-1)*self.F + f) for f in range(self.F)]
        self.Omega[usr] = [f for f in req_sub_files if f not in self.cache_contents[usr]]
        self.make_intervals()
        self.build_user_groups()

    def compute_eta(self, user_group, t):
        eta = 0
        for k in user_group:
            if (self.deadlines[k] == t):
                return float('inf')
            eta += (len(self.Omega[k]) - self.user_benefits[k]) / (self.deadlines[k] - t)
        return eta

    def run_online(self, all_arrival_times, all_deadlines, req_files, cache_contents):
        rate = 0
        self.user_benefits = [0 for _ in range(self.K)]
        inv_n = 1 / self.n
        for k in range(self.K):
            self.user_arrived(k, all_arrival_times[k], all_deadlines[k], req_files[k], cache_contents[k])
            LP_model = LpProblem('LP_DCAsynchronous', LpMinimize)
            #   LP variables
            x = LpVariable.matrix('x_var', list(range(self.num_var_x)), 0, None, LpContinuous)
            y = LpVariable.matrix('y_var', list(range(self.num_var_y)), 0, None, LpContinuous)
            #   First Constrain
            for nl in range(len(self.len_time_intervals)):
                tmp_const = self.constraint_one[nl]
                if len(tmp_const) > 0:
                    LP_model += lpSum([x[i] for i in tmp_const]) <= self.len_time_intervals[nl]
            #   second constraint
            for nug in range(len(self.label_usr_grp)):
                label_usr_grp = self.label_usr_grp[nug]
                sub_val = 0
                if label_usr_grp in self.submitted_usr_grp:
                    sub_val = self.submitted_usr_grp[label_usr_grp]
                for tmp_const in self.constraint_two_left[nug]:
                    if len(tmp_const) > 0:
                        LP_model += lpSum([y[j] for j in tmp_const]) <= lpSum([x[i] for i in self.constraint_two_right[nug]]) + sub_val
            #   third constraint
            for kk in self.active_users:
                for f in self.Omega[kk]:
                    LP_model += lpSum([y[i] for i in self.constraint_three[f]]) == 1

            #   cost function
            LP_model += lpSum(x), 'Min_Rate'

            #   solve LP
            LP_model.solve()
            if LpStatus[LP_model.status] == 'Infeasible':
                return None, -1

            #   interpert solution
            strt_time = all_arrival_times[k]
            if k < self.K-1:
                end_time = all_arrival_times[k+1]
            else:
                end_time = max(all_deadlines)
            opt_solutions = []

            for i in range(self.num_var_x):
                if value(x[i]) >= inv_n:
                    opt_solutions.append(XU(self.U_x[i], self.l_x[i], value(x[i])))

            # sort solution
            opt_solutions.sort()
            # pick solution
            t = strt_time
            while t < end_time:
            # for t in range(strt_time, end_time, inv_n):
                for ns in range(len(opt_solutions)):
                    tmp_sol = opt_solutions[ns]
                    if self.compute_eta(tmp_sol.user_group, t) >= self.eta_threshold:
                        label_usr_grp = ','.join(str(_) for _ in tmp_sol.user_group)
                        print('user group: %s at time %1.2f' % (label_usr_grp, t))
                        if label_usr_grp not in self.submitted_usr_grp:
                            self.submitted_usr_grp[label_usr_grp] = 0
                        opt_solutions[ns].val -= inv_n
                        if opt_solutions[ns].val < inv_n:
                            del opt_solutions[ns]
                        self.submitted_usr_grp[label_usr_grp] += inv_n
                        rate += inv_n
                        for i in tmp_sol.user_group:
                            self.user_benefits[i] += inv_n
                            if i in self.active_users and self.user_benefits[i] >= len(self.Omega[i]) - 1e-4:
                                self.active_users.remove(i)
                        break
                t += inv_n
        return self.submitted_usr_grp, rate