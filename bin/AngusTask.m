function tempVect = AngusTask(sMat)

sortedsMat = sort(sMat,2);
diffAdj = diff(sortedsMat,[],2);
diffAdj(diffAdj > 1) = 0;
% b = diffAdj == 1;

tempVect = zeros(size(sMat(:,1)))

for rowS = 1:size(diffAdj,1) %control statement
    if sum(diffAdj(rowS,:))>=1
        tempVect(rowS,1) = 1
    end
end


%I will need to use diff?
%At each row, get the 1 x m vector of that row.
%Within this vector, go through each element in order (starting from the
%second element.
%True and that kind of stuff?
%If none of the elements is exactly 1 greater
%than its previous (e.g., they are all
%2 greater than the previous), then
%fill the corresponding row in
%the output matrix with a 0.


% for participant = 1:length(Participants)
%     for Session = 1:length(Sessions)
%         Cactus = load(Participants{participant},Sessions{Session});
%         Values(participant,Session) = Cactus.(Sessions{Session}); %syntax for structures
%     end
% end


