% initialization tells the compiler what to execute as the initial goal.
% Here we are saying that the goal is to execute the main function.
% Basically this is what gets ran before the interpreter has a chance to do anything.
% http://www.gprolog.org/manual/gprolog.html#initialization%2F1
:- initialization(main).

% Straightforward stuff
fizzbuzz(X) :- 0 is X mod 15, write('FizzBuzz'), nl.
fizzbuzz(X) :- 0 is X mod 3, write('Fizz'), nl.
fizzbuzz(X) :- 0 is X mod 5, write('Buzz'), nl.
fizzbuzz(X) :- write(X), nl.

% This line says that when we get to doFizzBuzz(1) that we DO NOT backtrack.
% Delete this line to see what happens when backtracking occurs.
doFizzBuzz(0) :- !.

% Here we are recursively counting down from x, and performing FizzBuzz on X.
doFizzBuzz(X) :- 
    X2 is X - 1,
    fizzbuzz(X),
    doFizzBuzz(X2).

% The main function we defined in our initialization call.
main :- 
        
        % grabs the arguments passed in the command line as a list.
        % http://www.gprolog.org/manual/gprolog.html#hevea_default841
        argument_list(ARGS),

        % nth0 is a prolog list processing function to grab the nth element from a list.
        % http://www.gprolog.org/manual/html_node/gprolog044.html#sec221
        nth0(0, ARGS, Var),

        write(Var),
        
        % number_atom converts a "Number" to an atom.
        % http://www.gprolog.org/manual/html_node/gprolog043.html#number-atom%2F2
        number_atom(AsInt, Var),
        doFizzBuzz(AsInt),

        % makes it so that the interpreter immediately closes.
        % http://www.gprolog.org/manual/gprolog.html#abort%2F0
        halt(0).
