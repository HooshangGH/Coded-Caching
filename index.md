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

