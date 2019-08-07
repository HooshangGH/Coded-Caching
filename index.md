<head>
  <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.10.2/dist/katex.min.css" integrity="sha384-yFRtMMDnQtDRO8rLpMIKrtPCD5jdktao2TV19YiZYWMDkUR5GQZR/NOVTdquEx1j" crossorigin="anonymous">
<script defer src="https://cdn.jsdelivr.net/npm/katex@0.10.2/dist/katex.min.js" integrity="sha384-9Nhn55MVVN0/4OFx7EE5kpFBPsEMZxKTCnA+4fqDmg12eCTqGi6+BB2LjY8brQxJ" crossorigin="anonymous"></script>
<script defer src="https://cdn.jsdelivr.net/npm/katex@0.10.2/dist/contrib/auto-render.min.js" integrity="sha384-kWPLUVMOks5AQFrykwIup5lo0m3iMkkHrD0uJ4H5cjeGihAutqP0yW0J6dpFiVkI" crossorigin="anonymous" onload="renderMathInElement(document.body);"></script>
  </head>

# Asynchronous Coded Caching

## Original Coded Caching Problem
In their pioneering work, Maddah-Ali and Niesen considered the usage of coding in the caching problem. In this so-called "coded caching" setting, there is a server containing $N$ files, $K$ users each with a cache that can store up to $M$ files. The users are connected to the server via an error-free shared link.

The system operates in two distinct phases. 
* In the \textbf{placement phase} the content of the caches is populated by server. This phase does not depend on the future requests of the users which are assumed to be arbitrary. 
* In the \textbf{delivery phase} each user makes a request and the server transmits potentially coded signals to satisfy the requests of the users. The work of Maddah-Ali and Niesen demonstrated that significant reductions in the network traffic were possible as compared to conventional caching. Crucially, these gains continue to hold even if the popularity of the files is not taken into account.


## Asynchronous Coded Caching
While the coded caching scheme introduced by Maddah-Ali and Niesen achieves a significant result, the original formulation of the coded caching problem assumes that the user requests are synchronized, i.e., all file requests from the users arrive at the server at the same time. From a practical perspective, it is important to consider the case when the requests of the users are not synchronized; we refer to this as the asynchronous coded caching problem. 
<p align="center">
  <img src="Netflix.jpg" width="300" height="250">
</p>



