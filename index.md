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
* In the \textbf{placement phase} the content of the caches is populated by server. This phase does not depend on the future requests of the users which are assumed to be arbitrary. 
* In the \textbf{delivery phase} each user makes a request and the server transmits potentially coded signals to satisfy the requests of the users.

The work of Maddah-Ali and Niesen demonstrated that significant reductions in the network traffic were possible as compared to conventional caching. Crucially, these gains continue to hold even if the popularity of the files is not taken into account.


## Asynchronous Coded Caching
While the coded caching scheme introduced by Maddah-Ali and Niesen achieves a significant result, the original formulation of the coded caching problem assumes that the user requests are synchronized, i.e., all file requests from the users arrive at the server at the same time. From a practical perspective, it is important to consider the case when the requests of the users are not synchronized; we refer to this as the asynchronous coded caching problem. 
<p align="center">
  <img src="Netflix.jpg" width="300" height="250">
</p>



