//Solve   adapted from http://home.att.net/~srschmitt/
//to do: transpose, copyOnWrite
//saves data in row-major order (transpose of Matlab)
function Matrix(rows,cols,data)
{
 if (typeof rows == "string") // "0 1, 1 0"
 {
  var s = rows.trim().split(/[,\n]/g);
  this.rows = s.length;
  this.data = rows.replace(/(^\s+|\s+$)/g,'').split(/[,\s]+/g);
  this.cols = Math.floor(this.data.length / this.rows);
  for (var i = this.data.length - 1; i >= 0; i--) 
  {
    this.data[i] = Number(this.data[i]);
  }
 }
 else if (typeof rows === "object" && (rows.length || false))
 {
   this.rows = rows.length;
   this.cols = typeof rows[0] === "object" ? rows[0].length : 1;
   this.data = new Array(this.rows*this.cols);
    for (var r=this.rows-1; r>=0; r--)
      for (var c=this.cols-1; c>=0; c--)
        this.data[r*this.cols + c] = rows[r][c];
 }
 else
 {
   this.rows = rows;
   this.cols = cols;
   if (typeof data === "undefined") this.data = new Array(rows*cols);
   else this.data = data;
 }
 this.length = this.data.length;
}

Matrix.column = function(d)
{
  return new Matrix(d.length,1,d);
}

Matrix.row = function(d)
{
  return new Matrix(1,d.length,d);
}

Matrix.prototype.size = function()
{
  return this.rows * this.cols;
}

Matrix.prototype.column = function(c)
{
  var d = new Array(this.rows);
  for (var r=0; r<this.rows; r++) d[r] = this.data[r*this.cols + c];
  return new Matrix(this.rows,1,d);
}

Matrix.prototype.row = function(r)
{
  var i = r*this.cols;
  var d = this.data.slice(i,i+this.cols);
  return new Matrix(1,this.cols,d);
}

Matrix.prototype.toSource = function()
{
 return '(new Matrix('+this.rows+','+this.cols+','+this.data.toSource()+'))';
}

Matrix.prototype.zero = function()
{
 var max = this.rows *  this.cols;
 for (var r=0; r<max; r++)
  this.data[r] = 0;
return this;
}

Matrix.identity = function(s)
{
 m = new Matrix(s,s)
 m.zero();
 for (var r=s-1; r >=0; r--)
  m.set(r,r,1);
 return m;
}

Matrix.zeros = function(r,c)
{
 m = new Matrix(r,c || r)
 return m.zero();
}

Matrix.diagonal = function(s)
{
 m = new Matrix(s.length,s.length)
 m.zero();
 for (var r=0; r < m.rows; r++)
  m.set(r,r,s[r])
 return m;
}

Matrix.prototype.average = function(columnWise) //defaults row-wise average
{
  var n = columnWise ? this.cols: this.rows;
  var m = columnWise ? this.rows: this.cols;
  var data = new Array(m);
  for (var i = m - 1; i >= 0; i--) data[i]=0;
    for (var r=this.rows-1; r>=0; r--)
      for (var c=this.cols-1; c>=0; c--)
        data[columnWise ? r: c] += this.data[r*this.cols + c] / n;
  return new Matrix(columnWise? this.rows:1, columnWise?1:this.cols,data);
}

Matrix.prototype.add = function(m,s)
{
  var data;
  if (typeof s === 'undefined') s=1.0;
  if (m.constructor.name === 'Array' && m.length === this.data.length)
  {
    data = this.data.concat();
    for (var i=0; i<data.length; i++) data[i] += m[i];
    return data;
  }

  if (this.rows === m.rows && this.cols === m.cols)
  {
    data = new Array(this.rows*this.cols);
    for (var i = this.data.length - 1; i >= 0; i--) data[i] = this.data[i] + s*m.data[i];
  }
  else if (m.rows == this.rows && m.cols === 1) //m is a column, add column-wise
  {
    data = new Array(this.rows*this.cols);
    var i=data.length-1;
    for (var r=this.rows-1; r>=0; r--)
      for (var c=this.cols-1; c>=0; c--)
        data[i--] = this.data[r*this.cols + c] + s*m.data[r];
  }
  else if (m.cols == this.cols && m.rows === 1) //m is a row, add row-wise
  {
    data = new Array(this.rows*this.cols);
    var i=data.length-1;
    for (var r=this.rows-1; r>=0; r--)
      for (var c=this.cols-1; c>=0; c--)
        data[i--] = this.data[r*this.cols + c] + s*m.data[c];
  }
  else
  {
    throw "invalid matrix addition "+this.rows+'x'+this.cols+" "+m.rows+'x'+m.cols+"(try a.add(b) instead of b.add(a))"
  }
  return new Matrix(this.rows,this.cols,data);
}

Matrix.prototype.addTo = function(m)
{
  if (typeof m === 'number')
  {
        for (var i = this.data.length - 1; i >= 0; i--) this.data[i] += m;
  }
  if (this.rows === m.rows && this.cols === m.cols)
  {
    for (var i = this.data.length - 1; i >= 0; i--) this.data[i] += m.data[i];
  }
  else if (m.rows == this.rows && m.cols === 1) //m is a column, add column-wise
  {
    for (var r=this.rows-1; r>=0; r--)
      for (var c=this.cols-1; c>=0; c--)
        this.data[r*this.cols + c] += m.data[c];
  }
  else if (m.cols == this.cols && m.cols === 1) //m is a column, add column-wise
  {
    for (var r=this.rows-1; r>=0; r--)
      for (var c=this.cols-1; c>=0; c--)
        this.data[r*this.cols + c] += m.data[r];
  }
  return this;
}

//indexes start at zero
Matrix.prototype.at = function (r,c)
{
 return this.data[r*this.cols + c];
}

Matrix.prototype.get = function (r,c)
{
 return this.data[r*this.cols + c];
}

Matrix.prototype.put = function (r,c,x)
{
 this.data[r*this.cols + c]=x;
}

Matrix.prototype.set = function (r,c,x)
{
 this.data[r*this.cols + c]=x;
}

Matrix.prototype.toString = function(prec)
{
 var s = "";
 var x;
 for (var r=0; r<this.rows; r++)
 {
  for (var c =0; c<this.cols; c++)
 {
   x = this.at(r,c);
   s += (prec ? x.toFixed(prec) : x) + '\t'
 }
  s += '\n'
 }
 return s
}

Matrix.prototype.copy = function()
{
 return new Matrix(this.rows,this.cols,this.data.concat());
}

Matrix.prototype.transpose = function ()
{
 var data = new Array(this.rows*this.cols);
 var i=0;
 for (var c =0; c<this.cols; c++)
  for (var r=0; r<this.rows; r++)
  {
    data[i++] = this.data[r*this.cols + c];
 }
 return new Matrix(this.cols,this.rows,data);
}

Matrix.prototype.times = function (other)
{
 data = [];

 if ((typeof other) === "number") //multiplication
 {
   for (var i = 0; i < this.data.length; i++)
      data[i]=this.data[i] * other;
   return new Matrix(this.rows,this.cols,data);
 }

 if ((typeof other) === "object" && other.constructor !== this.constructor && other.length == this.cols)
 {
   for (var r=0; r<this.rows; r++)
   {
     var x = 0;
     for (var d =0; d<this.cols; d++)
       x += this.at(r,d) * other[d];
     data.push(x);
   }
   return data;
 }

 if (this.cols !== other.rows) 
  throw "Matrix mismatch " + this.size() + " " + other.size();
 for (var r=0; r<this.rows; r++)
  for (var c =0; c<other.cols; c++)
  {
   var x = 0;
   for (var d =0; d<this.cols; d++)
     x += this.at(r,d) * other.at(d,c)
   data.push(x);
 }

 return new Matrix(this.rows,other.cols,data);
}

Matrix.prototype.trace = function ()
{
 var r = 0;
 if (this.rows != this.cols) throw "Not square"

 for (var i = 0; i<this.rows; i++)
  r += this.at(i,i);
 return r;
}

Matrix.prototype.dot = function (other)
{
 var r = 0;
 if (this.data.length != 3 || other.data.length != 3) throw "Matrix mismatch"

 for (var i=0; i<this.data.length; i++)
   r += this.data[i] * other.data[i];
 return r;
}

Matrix.prototype.cross = function (other)
{
 var r = 0;
 if (this.data.length != 3 || other.data.length != 3) throw "Not 3-vector"

 return new Matrix(3,1,[this.data[1]*other.data[2]-this.data[2]*other.data[1],this.data[2]*other.data[0]-this.data[0]*other.data[2],this.data[0]*other.data[1]-this.data[1]*other.data[0]]);
}

Matrix.prototype.norm = function (other)
{
 var r = 0;
 if (this.rows != 1 && this.cols != 1) throw "Not a vector"
 for (var i=0; i<this.data.length; i++)
  r += this.data[i] * this.data[i];
 return Math.sqrt(r);
}

Matrix.prototype.det3 = function ()
{
 var r = 0;
 if (this.rows != 3 || this.cols != 3) throw "Not 3x3"

 return this.at(0,0)*(this.at(1,1)*this.at(2,2)-this.at(1,2)*this.at(2,1))-
    this.at(0,1)*(this.at(1,0)*this.at(2,2)-this.at(1,2)*this.at(2,0))-
    this.at(0,2)*(this.at(1,1)*this.at(2,0)-this.at(1,0)*this.at(2,1));
}

Matrix.prototype.scale = function (other)
{
 var data = []
 if (this.cols != other.cols || this.rows != other.rows) throw "Matrix mismatch";

 for (var i = 0; i<this.data.length; i++) data[i]=this.data[i] * other.data[i];
 return new Matrix(this.rows,this.cols,data);
}

Matrix.prototype.reciprocal = function ()
{
 var data = []
 for (var i = 0; i<this.data.length; i++) data[i] = (this.data[i] != 0) ? 1/this.data[i] : 0.0;
 return new Matrix(this.rows,this.cols,data);
}

Matrix.prototype.swapRows = function( row1, row2 )
{
 var t
 for (var i = 0; i< this.cols; i++)
 {
  t = this.get(row1, i)
  this.set(row1, i, this.get(row2, i))
  this.set(row2, i , t)
 }
}

Matrix.prototype.swapCols = function( col1, col2 )
{
 var t
 for (var i = 0; i< this.rows; i++)
 {
  t = this.get(i, col1)
  this.set(i, col1, this.get(i, col2))
  this.set(i, col2 , t)
 }
}

Matrix.prototype.round = function(precision)
{
  for (var i = 0; i<this.data.length; i++)
   data[i] = Math.round(this.data[i])
 return new Matrix(this.rows,this.cols,data);
}

Matrix.prototype.inverse = function(precision)
{//Gauss-Jordan elimination for a square matrix
   if (this.rows != this.cols) throw "Not square"

    var t = Matrix.identity(this.cols)
    var m = this.copy()
    var rank = this.rows

    function scale(m, row, value)
    {
      for (var i = 0; i< m.cols; i++)
        m.set(row, i, m.get(row, i) * value)
    }
    function combineRows(m, a, b, value)
    {
      for (var i = 0; i< m.cols; i++)
        m.set(a, i, m.get(a, i) + (value * m.get( b, i)))
    }

    if (!precision) precision = 1.0E-12
    else precision = Math.abs(precision)

    var i,j,d,pos
    for (pos =0; pos < rank; pos++)
    {
        d = m.get(pos, pos)

        if (Math.abs(d) <= precision) //done
        {
            m.set(pos,pos,0)
            break
        }

        //scale rows
        d = 1/d
        scale(m, pos, d)
        scale(t, pos, d)

        //subtract down
        for (j = pos+1; j < m.rows; j++)
        {
          d = m.get(j, pos)
          if (Math.abs(d) > precision)
          {
           combineRows(m, j, pos, -d)
           combineRows(t, j, pos, -d)
          }
        }
    }

    //backsubstitution on upper triangle
    for (var j = rank - 1; j > 0; j--)
     for (var i = 0; i < j; i++)
     {
      d = m.get(i, j)
      combineRows(m, i, j, -d)
      combineRows(t, i, j, -d)
     }
    return t
}

Matrix.rotation = function(axis,angle)
{
    if (axis.constructor === Matrix) axis = axis.data;
    var cq   = Math.cos(angle);
    var sq   = Math.sin(angle);
    var omcq = 1-cq;
    function sq(x) {return x*x;}
    return new Matrix(3,3,[ omcq*axis[0]*axis[0]+cq, omcq*axis[0]*axis[1]-sq*axis[2],  omcq*axis[0]*axis[2]+sq*axis[1],
        omcq*axis[1]*axis[0]+sq*axis[2],  omcq*axis[1]*axis[1]+cq,  omcq*axis[1]*axis[2]-sq*axis[0],
        omcq*axis[2]*axis[0]-sq*axis[1],  omcq*axis[2]*axis[1]+sq*axis[0],  omcq*axis[2]*axis[2]+cq]);
}

/** Given a square matrix S and vector Y, returns X such that S . X = Y */
Matrix.solve = function(S,Y)
{
  //adapted from http://home.att.net/~srschmitt/
  //eliminate
  var i, j, k;
  N = S.rows;
  if (N != S.cols) throw "S must be a square matrix";
  if (Y.constructor !== this.constructor && Y.length == N) Y = new Matrix(Y.length,1,Y);
  if (N != Y.rows) throw "Y must be a " + N + "x1 matrix";

  var A = new Matrix(S.rows, S.cols + Y.cols)

  for (i = 0; i < S.rows; i++)
  {
    A.put(i, S.cols, Y.at(i,0))
    for (j = 0; j < S.cols; j++)
      A.put(i, j, S.at(i,j))
  }

  for (i = 0; i < N; i++)
  {
    // find row with maximum in column i
    var max_row = i;
    for (j = i; j < N; j++)
    {
     if (Math.abs(A.at(j,i)) > Math.abs(A.at(max_row,i)))
       max_row = j;
    }

    // swap max row with row i of [A:y]
    for (k = i; k < N + 1; k++)
    {
     var tmp       = A.at(i,k);
     A.put(i,k, A.at(max_row,k));
     A.put(max_row, k, tmp);
    }

    // eliminate lower diagonal elements of [A]
    for (j = i + 1; j < N; j++)
    {
      for (k = N; k > i; k--)
      {
       if (A.at(i,i) == 0.0)
         return false;
       else
         A.put(j,k, A.at(j,k) - A.at(i,k)*A.at(j,i)/A.at(i,i));
      }
    }
  }

  //substitute
  var X = new Matrix(Y.rows, Y.cols); //cols should be 1
  X.zero();

  for (j = N - 1; j >= 0; j--)
  {
    var sum = 0.0;
    for (k = j+1; k < N; k++)
      sum += A.at(j,k)*X.at(k,0);

    X.put(j, 0, (A.at(j,N) - sum)/A.at(j,j));
  }
return X;
}

Matrix.fit = function(order,data)
{
 //data = [ [x1, y1], [x2, y2], [x3,y3] ]
 var N = data.length
 if (!order) throw "Matrix.fit missing parameter: order"
 var s = new Matrix(order+1,order+1)
 s.zero()

 var x = new Array(2*order +1); //[0,0,0,0,0]
 for (var i=0; i< x.length; i++)
  x[i] = 0;

 var xy = new Array(order+1); //[0,0,0]
 for (var i=0; i < xy.length; i++)
  xy[i] = 0;

 for (var d in data)
 {
  var a = data[d][0]
  var b = 1
  var c = data[d][1]
  for (var j=0; j< x.length; j++)
  {
   x[j] += b
   b = b * a
  }
  for (var j=0; j< xy.length; j++)
  {
   xy[order-j] += c
   c = c * a
  }
 }

 var y = new Matrix(xy.length,1,xy)

 for (var i=0; i<s.rows; i++)
 for (var j=0; j<s.cols; j++)
 {
  s.put(i, j, x[x.length - 1 - i - j ])
 }

 var ret = Matrix.solve(s,y).data
 ret.reverse()
 return ret;
}

/*
Translated from http://stitchpanorama.sourceforge.net/Python/svd.py

Here is the test case (first example) from Golub and Reinsch

a= new Matrix(8,5,[
    22.,10., 2.,  3., 7.,
    14., 7.,10.,  0., 8.,
    -1.,13.,-1.,-11., 3.,
    -3.,-2.,13., -2., 4.,
     9., 8., 1., -2., 4.,
     9., 1.,-7.,  5.,-1.,
     2.,-6., 6.,  5., 1.,
     4., 5., 0., -2., 2.] )

var [u,w,vt] = a.svd()
print(w)

[35.32704346531138,19.999999999999996,19.595917942265416,0,0]

// the correct answer is (eigenvalues sorted in descending order)

print([Math.sqrt(1248.),20.,Math.sqrt(384.),0.,0.])

[35.327043465311391, 20.0, 19.595917942265423, 0.0, 0.0]

*/


Matrix.prototype.pinv = function(precision)
{
 var prec= precision || 1.e-15   // assumes double prec
 var [u, w, v] = this.svd(prec)
 var tolerance = prec * Math.max(this.rows,this.cols) * Math.max.apply(null,w)

 wp = new Matrix(w.length,w.length)
 wp.zero()

 for (var i=0;i<w.length;i++)
  wp.set(i,i,Math.abs(w[i]) >= tolerance ? 1/w[i] : 0)

 return v.times(wp).times(u.transpose())
}

// Returns U, W, V, where (this= U * W * VT)
Matrix.prototype.svd= function(precision)
{
//Compute the thin SVD from G. H. Golub and C. Reinsch, Numer. Math. 14, 403-420 (1970)
  var prec= precision || 1.e-15 //Math.pow(2,-52) // assumes double prec
  var tolerance= 1.e-64/prec
  if (1.0+prec <= 1.0) throw "Need a bigger precision"
  if (!(tolerance > 0.0))  throw "Need a bigger tolerance"
  var itmax= 50
  var c=0
  var i=0
  var j=0
  var k=0
  var l=0

  var u= this.copy()
  var m= this.rows

  var n= this.cols

  if (m < n) throw "Need more rows than columns"

  var e = new Array(n)
  var q = new Array(n)
  for (i=0; i<n; i++) e[i] = q[i] = 0.0
  var v = new Matrix(n,n)
  v.zero()

  function pythag(a,b)
  {
    a = Math.abs(a)
    b = Math.abs(b)
    if (a > b)
      return a*Math.sqrt(1.0+(b*b/a/a))
    else if (b == 0.0)
      return a
    return b*Math.sqrt(1.0+(a*a/b/b))
  }

  //Householder's reduction to bidiagonal form

  var f= 0.0
  var g= 0.0
  var h= 0.0
  var x= 0.0
  var y= 0.0
  var z= 0.0
  var s= 0.0

  for (i=0; i < n; i++)
  {
    e[i]= g
    s= 0.0
    l= i+1
    for (j=i; j < m; j++)
      s += (u.get(j,i)*u.get(j,i))
    if (s <= tolerance)
      g= 0.0
    else
    {
      f= u.get(i,i)
      g= Math.sqrt(s)
      if (f >= 0.0) g= -g
      h= f*g-s
      u.set(i,i, f-g)
      for (j=l; j < n; j++)
      {
        s= 0.0
        for (k=i; k < m; k++)
          s += u.get(k,i)*u.get(k,j)
        f= s/h
        for (k=i; k < m; k++)
          u.set(k,j, u.get(k,j) + f*u.get(k,i))
      }
    }
    q[i]= g
    s= 0.0
    for (j=l; j < n; j++)
      s= s + u.get(i,j)*u.get(i,j)
    if (s <= tolerance)
      g= 0.0
    else
    {
      f= u.get(i,i+1)
      g= Math.sqrt(s)
      if (f >= 0.0) g= -g
      h= f*g - s
      u.set(i,i+1, f-g)
      for (j=l; j < n; j++) e[j]= u.get(i,j)/h
      for (j=l; j < m; j++)
      {
        s=0.0
        for (k=l; k < n; k++)
          s += (u.get(j,k)*u.get(i,k))
        for (k=l; k < n; k++)
          u.set(j,k, u.get(j,k)+(s*e[k]))
      }
    }
    y= Math.abs(q[i])+Math.abs(e[i])
    if (y>x)
      x=y
  }

  // accumulation of right hand gtransformations
  for (i=n-1; i != -1; i+= -1)
  {
    if (g != 0.0)
    {
      h= g*u.get(i,i+1)
      for (j=l; j < n; j++)
        v.set(j,i, u.get(i,j)/h)
      for (j=l; j < n; j++)
      {
        s=0.0
        for (k=l; k < n; k++)
          s += (u.get(i,k)*v.get(k,j))
        for (k=l; k < n; k++)
          v.set(k,j,v.get(k,j)+ (s*v.get(k,i)))
      }
    }
    for (j=l; j < n; j++)
    {
      v.set(i,j, 0.0)
      v.set(j,i, 0.0)
    }
    v.set(i,i, 1.0)
    g= e[i]
    l= i
  }

  // accumulation of left hand transformations
  for (i=n-1; i != -1; i+= -1)
  {
    l= i+1
    g= q[i]
    for (j=l; j < n; j++)
      u.set(i,j, 0.0)
    if (g != 0.0)
    {
      h= u.get(i,i)*g
      for (j=l; j < n; j++)
      {
        s=0.0
        for (k=l; k < m; k++) s += (u.get(k,i)*u.get(k,j))
        f= s/h
        for (k=i; k < m; k++) u.set(k,j,u.get(k,j) + (f*u.get(k,i)))
      }
      for (j=i; j < m; j++) u.set(j,i, u.get(j,i)/g)
    }
    else
      for (j=i; j < m; j++) u.set(j,i, 0.0)
    u.set(i,i,u.get(i,i) + 1.0)
  }

  // diagonalization of the bidiagonal form
  prec= prec*x
  for (k=n-1; k != -1; k+= -1)
  {
    for (var iteration=0; iteration < itmax; iteration++)
    { // test f splitting
      var test_convergence = false
      for (l=k; l != -1; l+= -1)
      {
        if (Math.abs(e[l]) <= prec)
        { test_convergence= true
          break
        }
        if (Math.abs(q[l-1]) <= prec)
          break
      }
      if (!test_convergence)
      { // cancellation of e[l] if l>0
        c= 0.0
        s= 1.0
        var l1= l-1
        for (i =l; i<k+1; i++)
        {
          f= s*e[i]
          e[i]= c*e[i]
          if (Math.abs(f) <= prec)
            break
          g= q[i]
          h= pythag(f,g)
          q[i]= h
          c= g/h
          s= -f/h
          for (j=0; j < m; j++)
          {
            y= u.get(j,l1)
            z= u.get(j,i)
            u.set(j,l1, y*c+(z*s))
            u.set(j,i, -y*s+(z*c))
          }
        }
      }
      // test f convergence
      z= q[k]
      if (l== k)
      { //convergence
        if (z<0.0)
        { //q[k] is made non-negative
          q[k]= -z
          for (j=0; j < n; j++)
            v.set(j,k, -v.get(j,k))
        }
        break  //break out of iteration loop and move on to next k value
      }
      if (iteration >= itmax-1)
        throw 'Error: no convergence.'
      // shift from bottom 2x2 minor
      x= q[l]
      y= q[k-1]
      g= e[k-1]
      h= e[k]
      f= ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
      g= pythag(f,1.0)
      if (f < 0.0)
        f= ((x-z)*(x+z)+h*(y/(f-g)-h))/x
      else
        f= ((x-z)*(x+z)+h*(y/(f+g)-h))/x
      // next QR transformation
      c= 1.0
      s= 1.0
      for (i=l+1; i< k+1; i++)
      {
        g= e[i]
        y= q[i]
        h= s*g
        g= c*g
        z= pythag(f,h)
        e[i-1]= z
        c= f/z
        s= h/z
        f= x*c+g*s
        g= -x*s+g*c
        h= y*s
        y= y*c
        for (j=0; j < n; j++)
        {
          x= v.get(j,i-1)
          z= v.get(j,i)
          v.set(j,i-1, x*c+z*s)
          v.set(j,i, -x*s+z*c)
        }
        z= pythag(f,h)
        q[i-1]= z
        c= f/z
        s= h/z
        f= c*g+s*y
        x= -s*g+c*y
        for (j=0; j < m; j++)
        {
          y= u.get(j,i-1)
          z= u.get(j,i)
          u.set(j,i-1, y*c+z*s)
          u.set(j,i, -y*s+z*c)
        }
      }
      e[l]= 0.0
      e[k]= f
      q[k]= x
    }
  }

  //vt= transpose(v)
  //return (u,q,vt)
  for (i=0;i<q.length; i++)
    if (q[i] < prec) q[i] = 0

  //sort eigenvalues
  for (i=0; i< n; i++)
  {
  //writeln(q)
   for (j=i-1; j >= 0; j--)
   {
    if (q[j] < q[i])
    {
  //  writeln(i,'-',j)
     c = q[j]
     q[j] = q[i]
     q[i] = c
     u.swapCols(i,j)
     v.swapCols(i,j)
     i = j
    }
   }
  }

  return [u,q,v];
}

Matrix.testpinv = function()
{
a= new Matrix(8,5,[
    22.,10., 2.,  3., 7.,
    14., 7.,10.,  0., 8.,
    -1.,13.,-1.,-11., 3.,
    -3.,-2.,13., -2., 4.,
     9., 8., 1., -2., 4.,
     9., 1.,-7.,  5.,-1.,
     2.,-6., 6.,  5., 1.,
     4., 5., 0., -2., 2.] )

//load('svd.js')
var [u,w,v]= a.svd()

writeln(w)
writeln([Math.sqrt(1248.),20.,Math.sqrt(384.),0.,0.])
writeln()

b = a.pinv()
writeln(b)
writeln(a.times(b))
//print(u,w,vt)
}

Matrix.test = function()
{
  i = Matrix.identity(3);
  writeln('i = Matrix.identity(3):');
  writeln(i);

  writeln('i.toSource():');
  writeln(i.toSource());
  writeln()

  writeln('eval(i.toSource()):');
  writeln(eval(i.toSource()));

  a = new Matrix(3,2,[1,0,1,1,1,2])
  b = new Matrix(3,1,[6,0,0])
  aT = a.transpose()
  writeln('a = new Matrix(3,2,[1,0,1,1,1,2])');
  writeln(a);

  writeln('b = new Matrix(3,1,[6,0,0])');
  writeln(b);

  writeln('Matrix.times(aT,b): 6 / 0 ')
  writeln(aT.times(b))

  writeln('Matrix.times(aT,a): 3 3 / 3 5')
  writeln(aT.times(a))

  writeln('Given a square matrix S and vector Y, find X such that S . X = Y')
  writeln('s = new Matrix(3,3,[1,1,-1,1,2,1,2,-1,1])')
  writeln('y = new Matrix(3,1,[0,8,3])')
  writeln('x = Matrix.solve(s,y)')
  s = new Matrix(3,3,[1,1,-1,1,2,1,2,-1,1])
  y = new Matrix(3,1,[0,8,3])
  x = Matrix.solve(s,y)
  writeln(x)

  writeln('Matrix.fit(): Curve fit to 4 x^2 + 50 x -40') // 1,0,1
  writeln('x = Matrix.fit(2,[[9,734],[3,146],[-12,-64]]')
  x = Matrix.fit(2,[[9,734],[3,146],[-12,-64]])
  writeln(x)

  writeln('Matrix.fit(): Curve fit to 2 x^3 + 5 x^2 + 7 x -13')
  x = Matrix.fit(3,[[5,397],[8,1387],[0,-13],[-12,-2833]]) //[[0,1],[-1,2],[1,2]])
  writeln(x)

  writeln()
  writeln('Inverse')
  s = new Matrix(3,3,[1,5,6,3,2,4,6,2,3])
  writeln(s)
  writeln(s.inverse())
  writeln(s.times(s.inverse()).round())
  writeln()

  s = new Matrix(3,3,[1,5,6,3,2,4,6,2,3])
  writeln(s)
  writeln(s.inverse())
  writeln(s.times(s.inverse()).round())
}

