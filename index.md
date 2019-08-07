<head>
  <script type="text/x-mathjax-config"> MathJax.Hub.Config({ TeX: { equationNumbers: { autoNumber: "all" } } }); </script>
       <script type="text/x-mathjax-config">
         MathJax.Hub.Config({
           tex2jax: {
             inlineMath: [ ['$','$'], ["\\(","\\)"] ],
             processEscapes: true
           }
         });
       </script>
       <script src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>

</head>

# Asynchronous Coded Caching

## Original Coded Caching Problem
In their pioneering work, Maddah-Ali and Niesen considered the usage of coding in the caching problem. In this so-called "coded caching" setting, there is a server containing $N$ files, $K$ users each with a cache that can store up to $M$ files. The users are connected to the server via an error-free shared link.

The system operates in two distinct phases. 
* In the *placement phase* the content of the caches is populated by server. This phase does not depend on the future requests of the users which are assumed to be arbitrary. 
* In the *delivery phase* each user makes a request and the server transmits potentially coded signals to satisfy the requests of the users.

The work of Maddah-Ali and Niesen demonstrated that significant reductions in the network traffic were possible as compared to conventional caching. Crucially, these gains continue to hold even if the popularity of the files is not taken into account.


## Asynchronous Coded Caching
While this is a significant result, the original formulation of the coded caching problem assumes that the user requests are synchronized, i.e., all file requests from the users arrive at the server at the same time. Henceforth, we refer to this as the synchronous setting. From a practical perspective, it is important to consider the asynchronous setting where user requests arrive at different times.
	In this case, a simple strategy would be to wait for the last request to arrive and then apply the scheme of Maddah-Ali and Niesen. Such a strategy will be quite good in terms of the overall rate of transmission from the server. However, this may be quite bad for an end user's experience, e.g., the delay experienced by the users will essentially be dominated by the arrival time of the last request.
	
<p align="center">
  <img src="Netflix.jpg" width="300" height="250">
</p>

In this work we formulate and study the coded caching problem when the user requests arrive at different times. Each user has a specific deadline by which his/her demand needs to be satisfied. The goal is to schedule transmission of packets so that each user is able to recover the requested file from the transmitted packets and his/her cache content within the prescribed deadline. We present algorithms for both the offline and online versions of this problem.

### Problem Formulation
We assume that time $\tau \geq 0$ is slotted. Let $[n]$ denote the set $\{1, \ldots, n\}$ and the symbol $\oplus$ represent the XOR operation. We assume that the server contains $N\geq K$ files\footnote{We assume that $N\geq K$ as it corresponds to the worst case rate where each of the $K$ users can request a different file. Furthermore, it is also the more practical scenario.} denoted by $W_{n}, n = 1, \dots, N$. The subfiles are denoted by $W_{n,f}$ so that $W_n = \{W_{n,f}: f \in [F]\}$ and the cache of user $i$ by $Z_i \subseteq \{ W_{n,f}: \ n \in [N], \ f \in [F] \}$. $Z_i$ contains at most $M_iF$ subfiles. In the delivery phase, user $i$ requests file $W_{d_i}$, where $d_i \in [N]$, from the server.
We let $\Omega^{(i)}$ denote the indices of the subfiles that are not present in the $i$-th user's cache, i.e.,
	\begin{align*}
	\Omega^{(i)} = \{ f:\ f \in [F], \ W_{d_i,f} \notin Z_i \}.
	\end{align*}
The equations in the delivery phase are assumed to be of the *all-but-one* type.
\begin{definition} *All-but-one equation*. Consider an equation $E$ such that
\begin{align*}
E = \oplus_{l=1}^\ell W_{d_{i_l}, f_{l}}.
\end{align*}
We say that $E$ is of the all-but-one type if for each $l \in [\ell]$, we have $W_{d_{i_l}, f_{l}} \notin Z_{i_l}$ and $W_{d_{i_l}, f_{l}} \in Z_{i_{k}}$ for all $k \in [\ell] \setminus \{l\}$.
\end{definition}
It is evident that an all-but-one equation transmitted from the server allows each of the users participating in the equation to recover a missing subfile that they need. The asynchronous coded caching problem can be formulated as follows.
	
*Inputs.*
* *User requests.* User $i$ requests file $W_{d_i}$, with $d_i \in [N]$ at time $T_i$. %User $i$'s request arrives at the server at time $T_i$.
* *Deadlines.* The $i$-th user needs to be satisfied by time $T_i + \Delta_i$, where $\Delta_i$ is a positive integer.
* *Transmission delay.* Each subfile needs $r$ time-slots to be transmitted over the shared link, i.e., each subfile can be treated as equivalent to $r$ packets, where each packet can be transmitted in one time slot.

As the problem is symmetric with respect to users, w.l.o.g. we assume that $T_1 \leq T_2 \leq \ldots \leq T_K$. Let $T_{\max} = \max_i (T_i + \Delta_i)$. Note that upon sorting the set of arrival times and deadlines, i.e., $\cup_{i=1}^{K} \{T_i, T_i + \Delta_i\}$, we can divide the interval $[T_1, T_{\max})$ into {\it at most} $2K-1$ non-overlapping intervals. Let the integer $\beta$, where $1 \leq \beta \leq 2K-1$ denote the number of intervals.
	
Let $\Pi_1, \ldots, \Pi_{\beta}$ represent the intervals where $\Pi_i$ appears before $\Pi_j$ if $i < j$; $|\Pi_\ell|$ denotes the length of interval $\Pi_\ell$ . The intervals are left-closed and right-open.  An easy to see but very useful property of the intervals that we have defined is that for a given $i$, either $[T_i, T_i + \Delta_i) \cap \Pi_\ell = \Pi_\ell$ or $[T_i, T_i + \Delta_i) \cap \Pi_\ell = \emptyset$. Fig. \ref{Fig:N_3_K_3_M_1_offline} shows an example when $K=3$. We define $U_\ell = \{i \in [K]:~ [T_i , T_i + \Delta_i) \cap \Pi_\ell = \Pi_\ell \}$, and $D_\ell = \{d_i\in [N]:~i\in U_\ell \}$.

Thus, $U_\ell$ is the set of active users in time interval $\Pi_\ell$ and $D_\ell$ is the corresponding set of active file requests.
	
*Outputs.*
* *Transmissions at each time slot.* If the problem is feasible, the schedule specifies which equations (of the all-but-one type) need to be transmitted at each time. The schedule is such that each user can recover all its missing subfiles within its deadline. The equations transmitted at time $\tau \in \Pi_\ell$ only depend on $D_\ell$. 


