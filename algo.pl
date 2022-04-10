:- initialization(main).
:- dynamic(input/1).
:- dynamic(fp/1).
:- dynamic(fm/1).
:- dynamic(tnh/1).
:- dynamic(mp/1).
:- dynamic(tns/1).
:- dynamic(outputfile/1).

%algorithm functions
    
allperms([],[]).

allperms(Nums, Result) :-
   setof(X, permutation(Nums, X), Result).

fmaggregator([],_,[]).
fmaggregator(Nums, [], Nums).
fmaggregator(Nums, [[I,J]|T], Res) :-
    nofm(Nums, I, J, X),
    fmaggregator(X, T, Res).
    
fpaggreator([],_,[]).
fpaggreator(Nums, [], Nums).
fpaggreator(Nums, [[I,J]|T], Res) :-
    yesfp(Nums, I, J, X),
    fpaggreator(X, T, Res).
    
tnhaggregator([],_,[]).
tnhaggregator(Nums, [], Nums).
tnhaggregator(Nums, [[I,J]|T], Res) :-
    notnh(Nums, I, J, X),
    tnhaggregator(X, T, Res).

    
%calculates all hard constraints
constraintcalc([],_,_,_,[]).
constraintcalc(Nums,Fp,Fm,Tnh,Res) :-
    fpaggreator(Nums, Fp, X),
    fmaggregator(X, Fm, Y),
    tnhaggregator(Y, Tnh, Res).
   
%performs check of entire list for 1 fm constraint   

nofm([],_,_,[]).

nofm([H|T], I, Num, Result) :-
    (   nth(I,H,Num)
    ->  Result = ResultT
    ;   Result = [H|ResultT]
    ),
    nofm(T,I,Num, ResultT).

%performs check of entire list for 1 fp constraint

yesfp([],_,_,[]).

yesfp([H|T], I, Num, Result) :-
    (   nth(I,H,Num)
    ->  Result = [H|ResultT]
    ;   Result = ResultT
    ),
    yesfp(T,I,Num, ResultT).

%performs check of entire list for 1 tnh constraint
    
notnh([],_,_,[]).

notnh([H|T], Num1, Num2, Result) :-
    (   tnhelper(Num1,Num2,H)
    ->  Result = ResultT
    ;   Result = [H|ResultT]
    ),
    notnh(T,Num1,Num2, ResultT).

tnhelper(Num1, Num2, [Num1,Num2|_]).
tnhelper(Num1, Num2, [Num2,_,_,_,_,_,_,Num1]). %change to [Num2,_,_,_,_,_,_,Num1].
tnhelper(Num1, Num2, [_|T]) :-
    tnhelper(Num1, Num2, T).

    
    
tnscalc(_,[],0).

tnscalc(Comb, [Num1,Num2,Val], Res) :-
    (   tnhelper(Num1,Num2,Comb)
    ->  Res = Val
    ;   Res = 0
    ).

    
    
mpcalc(_,[],0).

mpcalc([N1,N2,N3,N4,N5,N6,N7,N8], [One, Two, Three, Four, Five, Six, Seven, Eight], Res) :- %change to 1-8
    nth(N1, One, A),
    nth(N2, Two, B),
    nth(N3, Three, C),
    nth(N4, Four, D),
    nth(N5, Five, E),
    nth(N6, Six, F),
    nth(N7, Seven, G),
    nth(N8, Eight, H),
    Res is (A + B + C + D + E + F + G + H).

    
    
penaggregator(Comb,[],Mp,Val,Res) :-
    mpcalc(Comb, Mp, X),
    Res is (Val + X).

penaggregator(Comb, [Htns|T], Mp, Val, Res) :-
    tnscalc(Comb,Htns,X),
    Y is (Val + X),
    penaggregator(Comb, T, Mp, Y, Res).
    
    
 %prforms all penalty calculations and returns result. takes remainign combos x2, tns penalties, mp penalties, empty string, Result  
pencalc([],_,_,_,_,[[],[-1]]).

pencalc(Og,[H],Tns,Mp,Track,Result) :-
    penaggregator(H,Tns,Mp,0,X),
    append(Track,[X],Y),
    (   nth(1, Y, 0)
    ->  nth(1, Og, Combo), once(Result = [Combo, 0])
    ;   min_list(Y, Score), nth(N, Y, Score), nth(N, Og, Combo), once(Result = [Combo, Score])
    ).
    
pencalc(Og,[H|T],Tns,Mp,Track,Result) :-
    penaggregator(H,Tns,Mp,0,X),
    (   X = 0
    ->  Result = [H, 0]
    ;   append(Track,[X],Y), pencalc(Og,T,Tns,Mp,Y,Result)
    ).
    

    
    
%IO functions

writefile(File, Txt) :-
    open(File, write, Stream),
    write(Stream, Txt),
    close(Stream).
    
readfile(File) :-
    open(File, read, Stream),
    get_char(Stream, C1),
    processstream(C1, [], Stream),
    close(Stream).
    
processstream(end_of_file,Input,_) :-
    retractall(input(_)),
    assertz(input(Input)),
    !.
    
processstream(C, Input, Stream) :-
    append(Input, [C], Input2),
    get_char(Stream, C2),
    processstream(C2, Input2, Stream).

    
    
writesol(File, Combo, Sol) :-
  %  formatter(Comb, Fcomb), 
    open(File, write, Stream), 
    write(Stream, 'Solution '),
    write(Stream, Combo),
    write(Stream, '; Quality: '),
    write(Stream, Sol),
    close(Stream).

    
formatter([A,B,C,D,E,F,G,H], FormattedCombo) :-
   tochar(A, A1), tochar(B, B1), tochar(C, C1), tochar(D, D1), tochar(E, E1), tochar(F, F1), tochar(G, G1), tochar(H, H1),
   Sp = ' ',
   atom_concat(A1, Sp, As),
   atom_concat(As, B1, Ab),
   atom_concat(Ab, Sp, Bs),
   atom_concat(Bs, C1, Bc),
   atom_concat(Bc, Sp, Cs),
   atom_concat(Cs, D1, Cd),
   atom_concat(Cd, Sp, Ds),
   atom_concat(Ds, E1, De),
   atom_concat(De, Sp, Es),
   atom_concat(Es, F1, Ef),
   atom_concat(Ef, Sp, Fs),
   atom_concat(Fs, G1, Fg),
   atom_concat(Fg, Sp, Gs),
   atom_concat(Gs, H1, FormattedCombo).
   
tochar(1, 'A').
tochar(2, 'B').
tochar(3, 'C').
tochar(4, 'D').
tochar(5, 'E').
tochar(6, 'F').
tochar(7, 'G').
tochar(8, 'H').
 
%only use for tasks
toint('1', 1).
toint('2', 2).
toint('3', 3).
toint('4', 4).
toint('5', 5).
toint('6', 6).
toint('7', 7).
toint('8', 8).

%only use for machines
ismach('A', 1).
ismach('B', 2).
ismach('C', 3).
ismach('D', 4).
ismach('E', 5).
ismach('F', 6).
ismach('G', 7).
ismach('H', 8).

%only use for machine penalties
mpint('0', 0).
mpint('1', 1).
mpint('2', 2).
mpint('3', 3).
mpint('4', 4).
mpint('5', 5).
mpint('6', 6).
mpint('7', 7).
mpint('8', 8).
mpint('9', 9).

%Error handling functions
%Must be called in order: nowhitespace, checkheaders, namecheck, namenotnull, forcedpartialrem, forcedpartialdata,
%forbiddenmachinedata, toonearharddata, mpremover, (on og input: mpdatahelper, machinepenaltydata, bruh, bruh2, bruh3), toonearsoftdata
%-1 = no valid solution
%-2 = Error parsing input file
%-3 = invalid machine/task
%-4 = machine penalty error
%-5 = invalid penalty
%-6 = invalid task

nowhitespace(Input, Nospace) :-
    delete(Input, ' ', X),
    delete(X, '\n', Y),
    delete(Y, '\r', Z),
    delete(Z, '\r\n',Nospace).

removehead(T, 0, T).
    
removehead([_H|T], Num, Res) :-
    Y is (Num - 1),
    removehead(T, Y, Res).

remname([f,o,r,c,e,d,p,a,r,t,i,a,l,a,s,s,i,g,n,m,e,n,t,:|T],[f,o,r,c,e,d,p,a,r,t,i,a,l,a,s,s,i,g,n,m,e,n,t,:|T]).

remname([_H|T], Res) :-
    remname(T, Res).

checkheaders(Input, Res) :-
    (sublist(['N',a,m,e,:], Input), sublist([f,o,r,c,e,d,p,a,r,t,i,a,l,a,s,s,i,g,n,m,e,n,t,:], Input), sublist([f,o,r,b,i,d,d,e,n,m,a,c,h,i,n,e,:], Input), sublist([t,o,o,-,n,e,a,r,t,a,s,k,s,:], Input), sublist([m,a,c,h,i,n,e,p,e,n,a,l,t,i,e,s,:], Input), sublist([t,o,o,-,n,e,a,r,p,e,n,a,l,i,t,i,e,s], Input), Res = Input)
    ;   (sublist(['N',a,m,e,:], Input), sublist([f,o,r,c,e,d,p,a,r,t,i,a,l,a,s,s,i,g,n,m,e,n,t,:], Input), sublist([f,o,r,b,i,d,d,e,n,m,a,c,h,i,n,e,:], Input), sublist([t,o,o,-,n,e,a,r,t,a,s,k,s,:], Input), sublist([m,a,c,h,i,n,e,p,e,n,a,l,t,i,e,s,:], Input), Res = [-2])
    ;   (sublist(['N',a,m,e,:], Input), sublist([f,o,r,c,e,d,p,a,r,t,i,a,l,a,s,s,i,g,n,m,e,n,t,:], Input), sublist([f,o,r,b,i,d,d,e,n,m,a,c,h,i,n,e,:], Input), sublist([t,o,o,-,n,e,a,r,t,a,s,k,s,:], Input), Res = [-2])
    ;   (sublist(['N',a,m,e,:], Input), sublist([f,o,r,c,e,d,p,a,r,t,i,a,l,a,s,s,i,g,n,m,e,n,t,:], Input), sublist([f,o,r,b,i,d,d,e,n,m,a,c,h,i,n,e,:], Input), Res = [-2])
    ;   (sublist(['N',a,m,e,:], Input), sublist([f,o,r,c,e,d,p,a,r,t,i,a,l,a,s,s,i,g,n,m,e,n,t,:], Input), Res = [-2])
    ;   (sublist(['N',a,m,e,:], Input), Res = [-2])
    ;   Res = [-2].
    
    
namecheck(Input, Res) :-
    (   prefix(['N',a,m,e,:], Input)
    ->  removehead(Input, 5, Res)
    ;   Res = [-2]
    ).
    
namenotnull(Input, Res) :-
    (   prefix([f,o,r,c,e,d,p,a,r,t,i,a,l,a,s,s,i,g,n,m,e,n,t,:], Input)
    ->  Res = [-2]
    ;   remname(Input, Res)
    ).

forcedpartialrem(Input, Res) :-
    (   prefix([f,o,r,c,e,d,p,a,r,t,i,a,l,a,s,s,i,g,n,m,e,n,t,:], Input)
    ->  removehead(Input, 24, Res)
    ;   Res = [-2]
    ).

forcedpartialdata([f,o,r,b,i,d,d,e,n,m,a,c,h,i,n,e,:|T], Data, Res) :-
    retractall(fp(_)),
    assertz(fp(Data)),
    Res = T.

forcedpartialdata(['(',N,',',L,')'|T], Data, Res) :-
    (toint(N, X), ismach(L, Y), append(Data, [[X,Y]], Newdata), forcedpartialdata(T, Newdata, Res))
    ;   Res = [-3].

forcedpartialdata(_Wrong, _Data, [-2]).

forbiddenmachinedata([t,o,o,-,n,e,a,r,t,a,s,k,s,:|T], Data, Res) :-
    retractall(fm(_)),
    assertz(fm(Data)),
    Res = T.

forbiddenmachinedata(['(',N,',',L,')'|T], Data, Res) :-
    (toint(N, X), ismach(L, Y), append(Data, [[X,Y]], Newdata), forbiddenmachinedata(T, Newdata, Res))
    ;   Res = [-3].

forbiddenmachinedata(_Wrong, _Data, [-2]).

toonearharddata([m,a,c,h,i,n,e,p,e,n,a,l,t,i,e,s,:|T], Data, Res) :-
    retractall(tnh(_)),
    assertz(tnh(Data)),
    Res = T.

toonearharddata(['(',N,',',L,')'|T], Data, Res) :-
    (ismach(N, X), ismach(L, Y), append(Data, [[X,Y]], Newdata), toonearharddata(T, Newdata, Res))
    ;   Res = [-3].

toonearharddata(_Wrong, _Data, [-2]).

mpremover([t,o,o,-,n,e,a,r,p,e,n,a,l,i,t,i,e,s|T], T).

mpremover([_H|T], Res) :-
    mpremover(T, Res).
 
%must call this before machinepenaltydata, on original input with spaces
mpdatahelper([m,a,c,h,i,n,e,' ',p,e,n,a,l,t,i,e,s,:|T], Res) :-
    delete(T, '\r', X),
    mpnlremover(X, Res).

mpdatahelper([_H|T], Res) :-
    mpdatahelper(T, Res).
    
mpnlremover(['\n'|T], Res) :-
    mpnlremover(T, Res).
mpnlremover([' '|T], Res) :-
    mpnlremover(T, Res).
mpnlremover(Done, Done).

mpsizecheck(Data, 9, Data).
mpsizecheck(Data, I, Res) :-
    nth(I, Data, X),
    delete(X, [' '], Y),
    delete(Y, [], Z),
    delete(Z, ['\n'], A),
    (   length(A, 8)
    ->  J is (I + 1), mpsizecheck(Data, J, Res)
    ;   Res = [-4]
    ).

mpnumconverter(A, B) :-
    mpnumconverter(A, 0, B).    
mpnumconverter([], A, A).
mpnumconverter([H|T], A, N) :-
    X is (10*A + H),
    mpnumconverter(T, X, N).

mptonum([], Data, Data).    
mptonum([H|T], Data, Res) :-
    mptonumhelp(H, [], X),
    append(Data, [X], Newdata),
    mptonum(T, Newdata, Res).

mpcharint([], List, List).
mpcharint([H|T], List, Res) :-
    mpint(H, X),
    append(List, [X], Newlist),
    mpcharint(T, Newlist, Res).
    
mptonumhelp([], Row, Row).
mptonumhelp([H|T], Row, Res) :-
    mpcharint(H,[],I),
    mpnumconverter(I, X),
    append(Row, [X], Newrow),
    mptonumhelp(T, Newrow, Res).
    
mpisinthelp([], [0]).
mpisinthelp([H|T], Res) :-
    (   mpinter(H, X), X = [0]
    ->  mpisinthelp(T, Res)
    ;   Res = [-5]
    ).

mpinter([], [0]).
mpinter([_H|T], Res) :-
    (   mpint(_H, _X)
    ->  mpinter(T, Res)
    ;   Res = [-5]
    ).
    
mpisint([], [0]).
mpisint([H|T], Res) :-
    (   mpisinthelp(H, S), S = [-5]
    ->  Res = [-5]
    ;   mpisint(T, Res)
    ).

%must call these 3 in order after machinepenaltydata
bruh(Input, Res) :-
    (   length(Input, 8)
    ->  mpsizecheck(Input, 1, Res)
    ;   Res = [-4]
    ).    

bruh2(Input, Res) :-
    (   mpisint(Input, X), X = [0]
    ->  Res = Input
    ;   Res = [-5]
    ).
    
bruh3(Input, Res) :-
    mptonum(Input, [], X),
    retractall(mp(_)), 
    assertz(mp(X)),
    Res = X.
    
    
machinepenaltydata([t,o,o,-,n,e,a,r,' ',p,e,n,a,l,i,t,i,e,s|_T], _Thisnum, _Thisrow, Data, Res) :-
    delete(Data, [[]], Newdata),
    Res = Newdata.
%    (   length(Newdata, 8)
%    ->  mpsizecheck(Newdata, 1, X), %X = Newdata
%        ->  mpisint(Newdata, Y), Y = [-4]
%            ->  Res = [-65]
%            ;   mptonum(Newdata, [], Mp), retractall(mp(_)), assertz(mp(Mp)), Res = 0
 %       ;   Res = [-44]
%    ;   Res = [-34]
%    ).
    
machinepenaltydata(['\n'|T], Thisnum, Thisrow, Data, Res) :-
    append(Thisrow, [Thisnum], Newrow),
    append(Data, [Newrow], Newdata),
    machinepenaltydata(T, [], [], Newdata, Res).
    
machinepenaltydata([' '|T], Thisnum, Thisrow, Data, Res) :-
    append(Thisrow, [Thisnum], Newrow),
    machinepenaltydata(T, [], Newrow, Data, Res).

machinepenaltydata([N,'\n'|T], Thisnum, Thisrow, Data, Res) :-
    append(Thisnum, [N], Newnum),
    append(Thisrow, [Newnum], Newrow),
    append(Data, [Newrow], Newdata),
    machinepenaltydata(T, [], [], Newdata, Res).
    
machinepenaltydata([N, ' '|T], Thisnum, Thisrow, Data, Res) :-
    append(Thisnum, [N], Newnum),
    append(Thisrow, [Newnum], Newrow),
    machinepenaltydata(T, [], Newrow, Data, Res).   

machinepenaltydata([N,L|T], Thisnum, Thisrow, Data, Res) :-
    append(Thisnum, [N], Newnum),
    append(Newnum, [L], Newernum),
    machinepenaltydata(T, Newernum, Thisrow, Data, Res).
    
%machinepenaltydata(Wrong, Thisnum, Thisrow, Data, [-2]).
 
toonearsoftdata([], Data, Res) :-
    retractall(tns(_)),
    assertz(tns(Data)),
    Res = [].
    
toonearsoftdata(['(',N,',',L,',',P,')'|T], Data, Res) :-
    (ismach(N, X), ismach(L, Y), mpint(P, Z), append(Data, [[X,Y,Z]], Newdata), toonearsoftdata(T, Newdata, Res))
    ;   (ismach(N, X), ismach(L, Y), Res = [-5])
    ;   (ismach(N, X), Res = [-6])
    ;   Res = [-6].

toonearsoftdata(['(',N,',',L,',',P,Q,')'|T], Data, Res) :-
    (ismach(N, X), ismach(L, Y), mpint(P, Z1), mpint(Q, Z2), mpnumconverter([Z1,Z2], Z3), append(Data, [[X,Y,Z3]], Newdata), toonearsoftdata(T, Newdata, Res))
    ;   (ismach(N, X), ismach(L, Y), Res = [-5])
    ;   (ismach(N, X), Res = [-6])
    ;   Res = [-6].    

toonearsoftdata(['(',N,',',L,',',P,Q,R,')'|T], Data, Res) :-
    (ismach(N, X), ismach(L, Y), mpint(P, Z1), mpint(Q, Z2), mpint(R, Z3), mpnumconverter([Z1,Z2,Z3], Z4), append(Data, [[X,Y,Z4]], Newdata), toonearsoftdata(T, Newdata, Res))
    ;   (ismach(N, X), ismach(L, Y), Res = [-5])
    ;   (ismach(N, X), Res = [-6])
    ;   Res = [-6]. 
 
toonearsoftdata(_Wrong, _Data, [-2]).


%special error handlers
fpdoublebabysitter(Input, Res) :-
    (Input == [[]]; Res = [0])
    ;   fpdoubleshandler(Input, Res).

fpdoubleshandler([[_]|[]], [0]).
fpdoubleshandler([[]|_], [0]).
fpdoubleshandler([[I,J]|T], Res) :-
    (   fpdoubleshelper(I, J, T, X), X = [-7]
    ->  Res = [-7]
    ;   fpdoubleshandler(T, Res)
    ).

fpdoubleshelper(_I, _J, [], [0]).
fpdoubleshelper(I, J, [[I, J]|_T], [-7]).
fpdoubleshelper(I, _J, [[I, _K]|_T], [-7]).
fpdoubleshelper(_I, J, [[_K, J]|_T], [-7]).
fpdoubleshelper(I, J, [[_K, _L]|T], Res) :-
    fpdoubleshelper(I, J, T, Res).

%must be run through tns variables
ntsdupremover([],_,[]).
ntsdupremover([[I,J,K]|[]],Data, Res) :-
    append(Data, [[I,J,K]], X),
    retractall(tns(_)),
    assertz(tns(X)),
    Res = X.
    
ntsdupremover([[I,J,K]|T], Data, Res) :-
    (ntsduphelper(I,J,T,X), X = 0, append(Data, [[I,J,K]], Newdata), ntsdupremover(T, Newdata, Res))
    ;   (ntsdupremover(T, Data, Res)).
    
ntsduphelper(_I, _J, [], 0).
ntsduphelper(I, J, [[A, B, _K]|T], Res) :-
    (   I = A, J = B
    ->  Res = 1
    ;   ntsduphelper(I, J, T, Res)
    ).
 
%handlers for all error methods

fp2handler([-7]) :-
    outputfile(Out), writefile(Out, 'partial assignment error'), halt(0).
fp2handler(_I).

chhandler([-2]) :-
    outputfile(Out), writefile(Out, 'Error parsing input file1'), halt(0).
chhandler(_I).    
    
nchandler([-2]) :-
    outputfile(Out), writefile(Out, 'Error parsing input file2'), halt(0).
nchandler(_I).
    
nnnhandler([-2]) :-
    outputfile(Out), writefile(Out, 'Error parsing input file3'), halt(0).
nnnhandler(_I).

fprhandler([-2]) :-
    outputfile(Out), writefile(Out, 'Error parsing input file4'), halt(0).
fprhandler(_I).

fpdhandler([-3]) :-
    outputfile(Out), writefile(Out, 'invalid machine/task'), halt(0).
fpdhandler([-2]) :-
    outputfile(Out), writefile(Out, 'Error parsing input file5'), halt(0).
fpdhandler(_I).

fmdhandler([-3]) :-
    outputfile(Out), writefile(Out, 'invalid machine/task'), halt(0).
fmdhandler([-2]) :-
    outputfile(Out), writefile(Out, 'Error parsing input file6'), halt(0).
fmdhandler(_I).

tnhdhandler([-3]) :-
    outputfile(Out), writefile(Out, 'invalid machine/task'), halt(0).
tnhdhandler([-2]) :-
    outputfile(Out), writefile(Out, 'Error parsing input file7'), halt(0).
tnhdhandler(_I).
 
%dont use
mprhandler(Input, Res) :-
    (   mpremover(Input, X)
    ->  Res = X
    ;   outputfile(Out), writefile(Out, 'Error parsing input file8'), halt(0)
    ).
 
mpdhhandler(Input, Res) :-
    (   mpdatahelper(Input, X)
    ->  Res = X
    ;   outputfile(Out), writefile(Out, 'Error parsing input file9'), halt(0)
    ).


mpdhandler(Input, Res) :-
    (machinepenaltydata(Input, [], [], [], X)
    ->  Res = X
    ;   outputfile(Out), writefile(Out, 'Error parsing input file10'), halt(0)
    ).

b1handler([-4]) :-
    outputfile(Out), writefile(Out, 'machine penalty error'), halt(0).
b1handler(_I).

b2handler([-5]) :-
    outputfile(Out), writefile(Out, 'invalid penalty'), halt(0).
b2handler(_I).
    
%dont use    
b3handler(Input, Res) :-
    (   bruh3(Input, X)
    ->  Res = X
    ;   outputfile(Out), writefile(Out, 'Error parsing input file11'), halt(0)
    ).

tnsdhandler([-5]) :-
    outputfile(Out), writefile(Out, 'invalid penalty'), halt(0).
tnsdhandler([-6]) :-
    outputfile(Out), writefile(Out, 'invalid task'), halt(0).
tnsdhandler([-2]) :-
    outputfile(Out), writefile(Out, 'Error parsing input file12'), halt(0).
tnsdhandler(_I).
 
printhelp([_,[-1]]) :-
    outputfile(Out), writefile(Out, 'No valid solution possible!'), halt(0).
printhelp(_I).

resform([Combo|_T], Res) :-
    formatter(Combo, Res).
    
solform([_C,S], S).
    
 main :- 
    (argument_list(ARGS),
    nth(1, ARGS, Inputfile),
    nth(2, ARGS, Outputfile), retractall(outputfile(_)), assertz(outputfile(Outputfile)),
    readfile(Inputfile),
    input(A),
    nowhitespace(A, B),
    checkheaders(B, C),
    chhandler(C),
    namecheck(C, D),
    nchandler(D),
    namenotnull(D, E),
    nnnhandler(E),
    forcedpartialrem(E, F),
    fprhandler(F),
    forcedpartialdata(F, [], G),
    fpdhandler(G),
    forbiddenmachinedata(G, [], H),
    fmdhandler(H),
    toonearharddata(H, [], I),
    tnhdhandler(I),
    mpremover(I, J),
    mpdatahelper(A, Mp1),
    machinepenaltydata(Mp1, [], [], [], Mp2),
    bruh(Mp2, Mp3),
    b1handler(Mp3),
    bruh2(Mp3, Mp4),
    b2handler(Mp4),
    bruh3(Mp4, _Done),
    toonearsoftdata(J, [], K),
    tnsdhandler(K),     
    tns(L),
    ntsdupremover(L, [], _M),
    fp(N),
    fpdoublebabysitter(N, Doub),
    fp2handler(Doub),
    allperms([1,2,3,4,5,6,7,8], O),
    fm(Fm),
    fp(Fp),
    tnh(Tnh),
    mp(Mp),
    tns(Tns),
    constraintcalc(O, Fp, Fm, Tnh, P),
    pencalc(P, P, Tns, Mp, [], Result),
    printhelp(Result),
    resform(Result, Test),
    solform(Result, Test2),
    writesol(Outputfile, Test, Test2), halt(0))
    ;   write('Error'),
    halt(0).
    
    
    
    
    
    
