clc, clear all

% Lucca Rodrigues Pinto
% https://github.com/Lucca-Rodrigues-main

% Variaveis basicas
xb = [3 4];
% Variaveis nao basicas
xn = [1 2 5];
% Variaveis artificiais
fic = [5];
% Matriz de coeficientes A
A = [-1 -1 1 0 -1; 2 -1 0 1 -1];
% Coeficientes das basicas na funcao objetivo z
cb = [0 0];
% Coeficientes das basicas na funcao objetivo x0
cbL = [0 0];
% Coeficientes das nao basicas na funcao objetivo z
cn = [2 3 0];
% Coeficientes das nao basicas na funcao objetivo x0
cnL = [0 0 1];
% Coeficientes das basicas nas restricoes
B = A(:,xb);
% Coeficientes das nao basicas nas restricoes
N = A(:,xn);
% Vetor de resultados b
b = [-3; -2];

% Operacoes para construcao da tabela
op{1} = cb * inv(B) * N - cn;
op{2} = cbL * inv(B) * N - cnL;
op{3} = cb * inv(B) * b;
op{4} = cbL * inv(B) * b;
op{5} = inv(B) * N;
op{6} = inv(B) * b;
for i = 1:length(op)
    % Elimina imprecisao
    op{i}(abs(op{i}) < 1e-6) = 0;
end
fprintf('cb * inv(B) * N - cn =\n');
disp(op{1});
fprintf('cb'' * inv(B) * N - cn'' =\n');
disp(op{2});
fprintf('cb * inv(B) * b =\n');
disp(op{3});
fprintf('cb'' * inv(B) * b =\n');
disp(op{4});
fprintf('inv(B) * N =\n');
disp(op{5});
fprintf('inv(B) * b =\n');
disp(op{6});

% Gera a tabela inicial infactivel
[tab,xn,xb] = geratab2fasesAU(op, A, xb, xn, b);
tab

% Fase 1
while 1
    % Encontra o pivot
    [p,xn,xb] = pivotAUfase1(tab,xn,xb,fic);
    if all(p == 0)
        % Retornou pivot 0, fim da fase 1
        break
    end
    p
    disp('-------------------------');
    xn
    xb
    
    % Pivot vira 1
    tab(p(1),:) = tab(p(1),:)/tab(p(1),p(2));
    for i = 1:size(tab,1)
        if i ~= p(1)
            % Zerando demais valores na coluna do pivot
            tab(i,:) = tab(p(1),:)*(-tab(i,p(2))) + tab(i,:);
        end
    end
    tab
end
% Fase 2
% Elimina linha e coluna de x0
tab(2,:) = [];
tab(:,2) = [];
% Variaveis artificiais que nao estao na base
ficaux = fic;
ficaux(find(ismember(xb,fic))) = [];
% Elimina demais variaveis artificiais da tabela
tab(:,ficaux+1) = [];
xn(find(ismember(xn, ficaux))) = [];
disp('-------------------------');
tab
while 1
    % Encontra o pivot
    [p,xn,xb] = pivotAUfase2(tab,xn,xb,fic);
    if all(p == 0)
        % Retornou pivot 0, fim da fase 2
        break
    end
    p
    disp('-------------------------');
    xn
    xb
    % Pivot vira 1
    tab(p(1),:) = tab(p(1),:)/tab(p(1),p(2));
    for i = 1:size(tab,1)
        if i ~= p(1)
            % Zerando demais valores na coluna do pivot
            tab(i,:) = tab(p(1),:)*(-tab(i,p(2))) + tab(i,:);
        end
    end
    tab
end

function [tab,xn,xb] = geratab2fasesAU(op, A, xb, xn, b)
    % Gera a tabela inicial infactivel
    %    z  x0  xB        xN            RHS
    % z  1  0   0   cB B^-1 N - cN   cB B^-1 b
    % x0 0  1   0  cB' B^-1 N - cN'  cB' B^-1 b
    % xB 0  0   I       B^-1 N         B^-1 b
    
    % Inicializa com zero
    tab = sym(zeros(size(A,1)+2, size(A,2)+3));
    
    % Preenche as linhas de z e x0
    tab(1,1:2) = [1 0];
    tab(2,1:2) = [0 1];
    tab(1,[xn+2]) = op{1};
    tab(2,[xn+2]) = op{2};
    tab(1,end) = op{3};
    tab(2,end) = op{4};
    
    % Preenche as linhas de xB
    tab(3:size(A,1)+2,[xn+2]) = op{5};
    tab(3:size(A,1)+2,[xb+2]) = eye(length(xb));
    tab(3:end,end) = op{6};
    
    % Tornar a tabela factivel
    disp('-------------------------');
    disp('Tabela infactivel');
    
    % Pivotear a coluna de xa
    p(2) = size(tab,2)-1;
    % Valor mais negativo da coluna RHS linha xB
    [~,p(1)] = min(tab(3:end,end));
    p(1) = p(1) + 2;
    xn
    xb
    tab
    p
    
    disp('-------------------------');
    disp('Tabela factivel');
    % Pivot vira 1
    tab(p(1),:) = tab(p(1),:)/tab(p(1),p(2));
    for i = 1:size(tab,1)
        if i ~= p(1)
            % Zerando demais valores na coluna do pivot
            tab(i,:) = tab(p(1),:)*(-tab(i,p(2))) + tab(i,:);
        end
    end
    % Atualizando xB e xN
    iaux = xn(xn==p(2)-2);
    xn(xn==p(2)-2) = xb(p(1)-2);
    xb(p(1)-2) = iaux;
end

function [p,xn,xb] = pivotAUfase1(tab,xn,xb,fic)
    %    z  x0  xB        xN            RHS
    % z  1  0   0   cB B^-1 N - cN   cB B^-1 b
    % x0 0  1   0  cB' B^-1 N - cN'  cB' B^-1 b
    % xB 0  0   I       B^-1 N         B^-1 b
    
    if tab(2,end) == 0
        % Checa se RHS na linha de x0 atingiu valor 0
        if ~isempty(intersect(xb,fic)) && any(tab(find(ismember(xb,fic))+2,end) ~= 0)
            % Se houver variavel artificial na base com valor diferente de
            % zero, o problema pode nao ter solucao
            fprintf('problem may have no solution\n');
        else
            % Sem variaveis artificiais na base, fim da fase 1
            fprintf('end of phase 1\n');
        end
        p = [0 0];
        xn = xn;
        xb = xb;
        return
    else
        % Maximo valor positivo das nao basicas na linha de x0
        p(2) = max(tab(2,xn+2));
        % Armazena o indice
        iaux = find(tab(2,3:end-1)==p(2));
        iaux = iaux+2;
        if length(iaux) > 1
            % Checa se existe mais de um valor igual
            for i = 1:length(iaux)
                if ~isempty(find(tab(3:end,iaux(i)) > 0)) && all(ismember(fic,iaux(i)) == 0)
                    % Checa se existe ao menos um valor positivo na coluna
                    % da variavel nao basica selecionada e se
                    % preferivelmente nao eh uma variavel artificial
                    iaux = iaux(i);
                    break
                end
            end
            if length(iaux) > 1
                % Verifica se somente variaveis artificiais estao
                % disponiveis para entrar na base
                for i = 1:length(iaux)
                    if ~isempty(find(tab(3:end,iaux(i)) > 0))
                        % Checa se existe ao menos um valor positivo na
                        % coluna da variavel nao basica selecionada
                        iaux = iaux(i);
                        break
                    end
                end
                if length(iaux) > 1
                    % Nao existem valores positivos na coluna da variavel
                    % nao basica selecionada, entao o problema pode ser
                    % ilimitado
                    fprintf('Problem may be not limited\n');
                    p = [0 0];
                    xn = xn;
                    xb = xb;
                    return
                end
            end
        end
        if length(iaux) > 1 && ~isempty(intersect(iaux,fic))
            % Double checa se existe mais de um valor igual e se sao
            % variaveis artificiais
            iaux = intersect(iaux,fic);
        end
        p(2) = iaux(1);
        % Indices dos valores > 0 na coluna selecionada
        index = find(tab(3:end,p(2)) > 0);
        if isempty(index)
            % Checa se nao existem valores positivos na coluna selecionada
            fprintf('Problem may be not limited\n');
            p = [0 0];
            xn = xn;
            xb = xb;
            return
        end
        index = index + 2;
        % Indices dos menores valores de RHS dividido pelos valores
        % positivos da coluna selecionada
        iaux = find((tab(index,end)./tab(index,p(2))) == ...
            min(tab(index,end)./tab(index,p(2))));
        if length(iaux) > 1
            % Checa se existe mais de um valor igual
            for i = 1:length(iaux)
                if any(ismember(fic,iaux(i)) == 1)
                    % Se houver uma variavel artificial candidata a sair da
                    % base, sera selecionada
                    iaux = iaux(i);
                    break
                end
            end
        end
        p(1) = index(iaux(1));
        
        % Atualizando xB e xN
        iaux = xn(xn==p(2)-2);
        xn(xn==p(2)-2) = xb(p(1)-2);
        xb(p(1)-2) = iaux;
    end
end

function [p,xn,xb] = pivotAUfase2(tab,xn,xb,fic)
    %    z  xB        xN            RHS
    % z  1  0   cB B^-1 N - cN   cB B^-1 b
    % xB 0  I       B^-1 N         B^-1 b
    
    if all(tab(1,xn+1) <= 1e-6)
        % A fase 2 acaba quando os valores das variaveis nao basicas na
        % linha de z sao <= 0
        fprintf('end of phase 2\n');
        if any(tab(2:end,end) < 0)
            fprintf('primal problem is infeasible\n');
        end
        p = [0 0];
        xn = xn;
        xb = xb;
        return
    else
        % Maximo valor positivo das nao basicas na linha de z
        p(2) = max(tab(1,xn+1));
        % Armazena o indice
        iaux = find(tab(1,2:end-1)==p(2));
        iaux = iaux+1;
        if length(iaux) > 1
            % Checa se existe mais de um valor igual
            for i = 1:length(iaux)
                if ~isempty(find(tab(2:end,iaux(i)) > 0)) && all(ismember(fic,iaux(i)) == 0)
                    % Checa se existe ao menos um valor positivo na coluna
                    % da variavel nao basica selecionada e se
                    % preferivelmente nao eh uma variavel artificial
                    iaux = iaux(i);
                    break
                end
            end
            if length(iaux) > 1
                % Verifica se somente variaveis artificiais estao
                % disponiveis para entrar na base
                for i = 1:length(iaux)
                    if ~isempty(find(tab(2:end,iaux(i)) > 0))
                        % Checa se existe ao menos um valor positivo na
                        % coluna da variavel nao basica selecionada
                        iaux = iaux(i);
                        break
                    end
                end
                if length(iaux) > 1
                    % Nao existem valores positivos na coluna da variavel
                    % nao basica selecionada, entao o problema pode ser
                    % ilimitado
                    fprintf('Problem may be not limited\n');
                    p = [0 0];
                    xn = xn;
                    xb = xb;
                    return
                end
            end
        end
        p(2) = iaux(1);
        % Indices dos valores > 0 na coluna selecionada
        index = find(tab(2:end,p(2)) > 0);
        if isempty(index)
            % Checa se nao existem valores positivos na coluna selecionada
            fprintf('Problem may be not limited\n');
            p = [0 0];
            xn = xn;
            xb = xb;
            return
        end
        index = index + 1;
        % Indices dos menores valores de RHS dividido pelos valores
        % positivos da coluna selecionada
        iaux = find((tab(index,end)./tab(index,p(2))) == ...
            min(tab(index,end)./tab(index,p(2))));
        if length(iaux) > 1
            % Checa se existe mais de um valor igual
            for i = 1:length(iaux)
                if any(ismember(fic,iaux(i)) == 1)
                    % Se houver uma variavel artificial candidata a sair da
                    % base, sera selecionada
                    iaux = iaux(i);
                    break
                end
            end
        end
        p(1) = index(iaux(1));
        
        % Atualizando xB e xN
        iaux = xn(xn==p(2)-1);
        xn(xn==p(2)-1) = xb(p(1)-1);
        xb(p(1)-1) = iaux;
    end
end