designPositions=[  35 0 0; -17.5 30.3 0;-17.5 -30.3 0; 0 35 -70]';
offsetA = mean(designPositions,2);
disp(designPositions')
rotm = rotationMatrix(normr([1,2,3]),pi/3);
offset = [124,567,890]';
measuredPositions=rotm*designPositions + offset;
disp(measuredPositions')
N = size(designPositions,2);
weights = ones(1,N);
weights = 1/length(weights);

frameA = designPositions - offsetA;
offsetB = mean(measuredPositions,2);
frameB = measuredPositions - offsetB;

B = zeros(3,3);
for i=1:N % [3x1] * [1x3] gives the [3x3] outer product
  B = B + frameB(:,i) * frameA(:,i)';
end

[U,S,V] = svd(B);

M = diag([1 1 det(U)*det(V)]);

R = U * M * V';
disp('Inferred rotation')
disp(R)
angle = acos((trace(R)-1)/2)
% great, now how do I find the eigenvector for eivenvalue 1?

a = R-eye(3); 
b = a * [0 0 1]'
a(1:2,1:2)\b(1:2)
inferredPositions = R * (frameA) + offsetB;

%%
measuredDifference = measuredPositions - inferredPositions;

frameOffset = mean(measuredDifference,2);
disp('Coordinate system origin (in the measurement frame)')
disp(sprintf('%4.4g ',frameOffset))

pointErrors = measuredDifference - repmat(frameOffset,1,N);
disp('Point error vectors')
disp(pointErrors)

disp('Differential measurement errors')
disp(pointErrors - repmat(mean(pointErrors,2),1,N))

disp('Coordinate frame transformation error * 1000')
disp(1000*(R - rotm))

metric = norm( frameOffset);