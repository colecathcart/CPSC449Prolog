
%algorithm functions

perm([],[]).

perm(L,[H|T]) :-
    append(V,[H|U],L),
    append(V,U,W), 
    perm(W,T).


    
allperms([],[]).

allperms(Nums, Result) :-
   setof(X, perm(Nums, X), Result).

   

nofm([],_,_,[]).

nofm([H|T], I, Num, Result) :-
    (   nth(I,H,Num)
    ->  Result = ResultT
    ;   Result = [H|ResultT]
    ),
    nofm(T,I,Num, ResultT).



yesfp([],_,_,[]).

yesfp([H|T], I, Num, Result) :-
    (   nth(I,H,Num)
    ->  Result = [H|ResultT]
    ;   Result = ResultT
    ),
    yesfp(T,I,Num, ResultT).

    
    
notnh([],_,_,[]).

notnh([H|T], Num1, Num2, Result) :-
    (   tnhelper(Num1,Num2,H)
    ->  Result = ResultT
    ;   Result = [H|ResultT]
    ),
    notnh(T,Num1,Num2, ResultT).

tnhelper(Num1, Num2, [Num1,Num2|_]).
tnhelper(Num1, Num2, [Num2,_,Num1]). %change to [Num2,_,_,_,_,_,_,Num1].
tnhelper(Num1, Num2, [_|T]) :-
    tnhelper(Num1, Num2, T).

    
    
tnscalc(_,[],0).

tnscalc(Comb, [Num1,Num2,Val], Res) :-
    (   tnhelper(Num1,Num2,Comb)
    ->  Res = Val
    ;   Res = 0
    ).

    
    
mpcalc(_,[],0).

mpcalc([N1,N2,N3], [One, Two, Three], Res) :- %change to 1-8
    nth(N1, One, X),
    nth(N2, Two, Y),
    nth(N3, Three, Z),
    Res is (X + Y + Z).

    
    
penaggregator(Comb,[],Mp,Val,Res) :-
    mpcalc(Comb, Mp, X),
    Res is (Val + X).

penaggregator(Comb, [Htns|T], Mp, Val, Res) :-
    tnscalc(Comb,Htns,X),
    Y is (Val + X),
    penaggregator(Comb, T, Mp, Y, Res).
    
    
    
pencalc([],_,_,_,_,[[],[-1]]).

pencalc(Og,[H],Tns,Mp,Track,Result) :-
    penaggregator(H,Tns,Mp,0,X),
    append(Track,[X],Y),
    min_list(Y, Score),
    nth(N,Y,Score),
    nth(N,Og,Combo),
    once(Result = [Combo, Score]).
    
pencalc(Og,[H|T],Tns,Mp,Track,Result) :-
    penaggregator(H,Tns,Mp,0,X),
    append(Track,[X],Y),
    pencalc(Og,T,Tns,Mp,Y,Result).
    

    
    
%IO functions

writefile(File, Txt) :-
    open(File, write, Stream),
    write(Stream, Txt), nl,
    close(Stream).
    
readfile(File) :-
    open(File, read, Stream),
    get_char(Stream, C1),
    processstream(C1, Stream),
    close(Stream).
    
processstream(end_of_file,_) :- !.

processstream(C, Stream) :-
    write(C),
    get_char(Stream, C2),
    processstream(C2, Stream).
