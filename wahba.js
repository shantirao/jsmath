load('matrix.js')
load('geometry.js')


function getCorrdinateTransformation(points,apriori,measured)
{

	function listify(t)
	{
		var r;
	    if (t.indexOf(',') != -1)
		  r = t.trim().split(',');
		else r = t.trim().split(/[\s]+/g);
		if (r[r.length-1] == null) r.pop();
		return r.map(Number);
	}

// column vectors = XYZ
var designPositions = new Matrix(points.map(function(x){return x[apriori]})).transpose();
var measuredPositions = new Matrix(points.map(function(x){return x[measured]})).transpose();

 if (measuredPositions.length < 3 || measuredPositions.length != designPositions.length) return null;

 var offsetA = designPositions.average(true); //column vector
 var offsetB = measuredPositions.average(true);

 var frameA = designPositions.add(offsetA,-1);
 var frameB = measuredPositions.add(offsetB,-1);

 var B = Matrix.zeros(3);
 B.zero();

 for (var i=0; i<points.length; i++)
 {
  //  B.addTo(Matrix.column(measuredPositions[i]).times(Matrix.row(designPositions[i])));
  var x = frameB.column(i);
  var y = frameA.column(i).transpose();
  B.addTo(x.times(y));
  //row/col opposite of MATLAB because these are row-major
 }

 var [U,S,V] = B.svd();

//writeln('offsetA')
//writeln(offsetA)
//writeln('offsetB')
//writeln(offsetB)

 var M = Matrix.diagonal([1, 1, U.det3()*V.det3()]);

 var R = U.times(M).times(V.transpose());

 var prediction = R.times(frameA).add(offsetB);


 //var offset = offsetB.add(R.times(offsetA),-1);

 var differences = measuredPositions.add(prediction,-1);

 var euler = rotationToEuler(R); //
 writeln('Euler angles'); writeln(euler)
 writeln('Angles between axes')
writeln([0,1,2].map(function(x){return (180/Math.PI) * Math.asin(R.at(x,x))}));

//R is rotation matrix in the lab frame.
var z = new Matrix(3,1,[0,0,-1]);
var a = R.add(Matrix.identity(3),-1);
var Rz = a.times(z); 

var XY = new Matrix(2,2,[a.at(0,0),a.at(0,1),a.at(1,0),a.at(1,1)])
var axis = Matrix.solve(XY,new Matrix(2,1,Rz.data));
axis.data.push(1)
axis.rows++;
 

 return {differences: differences, offset: offsetB.add(R.times(offsetA),-1), rotation:R, euler: euler, angle:Math.acos((R.trace()-1)/2), axis:axis}
}