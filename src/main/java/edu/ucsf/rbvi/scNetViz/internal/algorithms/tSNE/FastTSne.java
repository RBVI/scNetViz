package edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE;

import static edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE.EjmlOps.addRowVector;
import static edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE.EjmlOps.assignAllLessThan;
import static edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE.EjmlOps.assignAtIndex;
import static edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE.EjmlOps.biggerThan;
import static edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE.EjmlOps.colMean;
import static edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE.EjmlOps.extractDoubleArray;
import static edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE.EjmlOps.maximize;
import static edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE.EjmlOps.replaceNaN;
import static edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE.EjmlOps.setData;
import static edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE.EjmlOps.setDiag;
import static edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE.EjmlOps.tile;
import static edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE.MatrixOps.abs;
import static edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE.MatrixOps.addColumnVector;
import static edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE.MatrixOps.assignValuesToRow;
import static edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE.MatrixOps.concatenate;
import static edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE.MatrixOps.equal;
import static edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE.MatrixOps.fillMatrix;
import static edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE.MatrixOps.getValuesFromRow;
import static edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE.MatrixOps.mean;
import static edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE.MatrixOps.negate;
import static edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE.MatrixOps.range;
import static edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE.MatrixOps.rnorm;
import static edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE.MatrixOps.scalarInverse;
import static edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE.MatrixOps.scalarMult;
import static edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE.MatrixOps.sqrt;
import static edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE.MatrixOps.square;
import static edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE.MatrixOps.sum;
import static edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE.MatrixOps.times;
import static org.ejml.dense.row.CommonOps_DDRM.add;
import static org.ejml.dense.row.CommonOps_DDRM.addEquals;
import static org.ejml.dense.row.CommonOps_DDRM.divide;
import static org.ejml.dense.row.CommonOps_DDRM.elementDiv;
import static org.ejml.dense.row.CommonOps_DDRM.elementExp;
import static org.ejml.dense.row.CommonOps_DDRM.elementLog;
import static org.ejml.dense.row.CommonOps_DDRM.elementMult;
import static org.ejml.dense.row.CommonOps_DDRM.elementPower;
import static org.ejml.dense.row.CommonOps_DDRM.elementSum;
import static org.ejml.dense.row.CommonOps_DDRM.mult;
import static org.ejml.dense.row.CommonOps_DDRM.multAddTransB;
import static org.ejml.dense.row.CommonOps_DDRM.scale;
import static org.ejml.dense.row.CommonOps_DDRM.subtract;
import static org.ejml.dense.row.CommonOps_DDRM.subtractEquals;
import static org.ejml.dense.row.CommonOps_DDRM.sumRows;
import static org.ejml.dense.row.CommonOps_DDRM.transpose;

import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;

import java.math.BigDecimal;

import org.ejml.data.DMatrixRMaj;

import org.cytoscape.work.TaskMonitor;
/**
*
* Author: Leif Jonsson (leif.jonsson@gmail.com)
* 
* This is a Java implementation of van der Maaten and Hintons t-sne 
* dimensionality reduction technique that is particularly well suited 
* for the visualization of high-dimensional datasets
*
*/
public class FastTSne implements TSne {
	MatrixOps mo = new MatrixOps();
	protected volatile boolean abort = false;
	static int DIGITS = 13;
	TaskMonitor monitor;
	TSneConfiguration config;

	public static double[][] readBinaryDoubleMatrix(int rows, int columns, String fn) throws FileNotFoundException, IOException {
		File matrixFile = new File(fn);
		double [][] matrix = new double[rows][columns];
		try (DataInputStream dis =
				new DataInputStream(new BufferedInputStream(new FileInputStream(matrixFile.getAbsolutePath())))) {
			for (int i = 0; i < matrix.length; i++) {
				for (int j = 0; j < matrix[0].length; j++) {
					matrix[i][j] = dis.readDouble();
				}
			}
		}
		return matrix;
	}
	
	//
	@Override
	public double [][] tsne(TSneConfiguration config, TaskMonitor monitor) {
		double[][] X      = config.getXin();
		int no_dims       = config.getOutputDims();
		int initial_dims  = config.getInitialDims(); 
		double perplexity = config.getPerplexity();
		int max_iter      = config.getMaxIter();
		boolean use_pca   = config.usePca();
		this.monitor      = monitor;
		this.config       = config;

		String IMPLEMENTATION_NAME = this.getClass().getSimpleName();
		// System.out.println("X:Shape is = " + X.length + " x " + X[0].length);
		monitor.showMessage(TaskMonitor.Level.INFO, "Running " + IMPLEMENTATION_NAME + ".");
		long end = System.currentTimeMillis();
		long start = System.currentTimeMillis();
		long total = System.currentTimeMillis();

		// Scale the data if we're supposed to
		if (config.logNormalize()) {
			X = MatrixOps.log(X, true);
		}

		if (config.centerAndScale()) {
			X = MatrixOps.centerAndScale(X);
		}
		System.gc();

		// Initialize variables
		if(use_pca && X[0].length > initial_dims && initial_dims > 0) {
			monitor.showMessage(TaskMonitor.Level.INFO, "Using PCA to reduce dimensions");
			PrincipalComponentAnalysis pca = new PrincipalComponentAnalysis();
			X = pca.pca(X, initial_dims);
			monitor.showMessage(TaskMonitor.Level.INFO, "X:Shape after PCA is = " + X.length + " x " + X[0].length);

		}

		// Since the original X is done, we should be able to
		// dispose of some memory
		System.gc();

		int n = X.length;
		double momentum = .5;
		double initial_momentum = 0.5;
		double final_momentum   = 0.8;
		int eta                 = 500;
		double min_gain         = 0.01;
		DMatrixRMaj Y        = new DMatrixRMaj(rnorm(n,no_dims));
		// DMatrixRMaj Y        = new DMatrixRMaj(fillMatrix(n,no_dims,.5));
		DMatrixRMaj Ysqlmul  = new DMatrixRMaj(Y.numRows,Y.numRows);
		DMatrixRMaj dY       = new DMatrixRMaj(fillMatrix(n,no_dims,0.0));
		DMatrixRMaj iY       = new DMatrixRMaj(fillMatrix(n,no_dims,0.0));
		DMatrixRMaj gains    = new DMatrixRMaj(fillMatrix(n,no_dims,1.0));
		DMatrixRMaj btNeg    = new DMatrixRMaj(n,no_dims);
		DMatrixRMaj bt       = new DMatrixRMaj(n,no_dims);

		if (config.cancelled())
			return null;

		writeMatrix("X", X);
		
		// Compute P-values
		DMatrixRMaj P        = new DMatrixRMaj(x2p(X, 1e-5, perplexity).P); // P = n x n
		writeMatrix("P1", P);

		// OK, now free up X
		X = null;
		System.gc();

		if (config.cancelled())
			return null;

		DMatrixRMaj Ptr      = new DMatrixRMaj(P.numRows,P.numCols);
		DMatrixRMaj L        = new DMatrixRMaj(P); // L = n x n
		DMatrixRMaj logdivide = new DMatrixRMaj(P.numRows,P.numCols);
		DMatrixRMaj diag     = new DMatrixRMaj(fillMatrix(L.numRows,L.numCols,0.0));
		
		transpose(P,Ptr);
		addEquals(P,Ptr);
		writeMatrix("P2", P);
		divide(P ,round(elementSum(P),DIGITS));
		replaceNaN(P,Double.MIN_VALUE);
		scale(4.0,P);					// early exaggeration
		maximize(P, 1e-12);

		if (config.cancelled())
			return null;
		
		monitor.showMessage(TaskMonitor.Level.INFO, "Y:Shape is = " + Y.getNumRows() + " x " + Y.getNumCols());

		DMatrixRMaj sqed  = new DMatrixRMaj(Y.numRows,Y.numCols);
		DMatrixRMaj sum_Y = new DMatrixRMaj(1,Y.numRows);
		DMatrixRMaj num   = new DMatrixRMaj(Y.numRows, Y.numRows);
		DMatrixRMaj Q     = new DMatrixRMaj(P.numRows,P.numCols);

		//System.out.println("Created sum_Y: ("+
		//                   ""+sum_Y.numRows+","+sum_Y.numCols+")");
		
		double progress = 0.0;
		for (int iter = 0; iter < max_iter && !abort; iter++) {
			if (config.cancelled())
				return null;
			progress = (double)iter/(double)max_iter;
			monitor.setProgress(progress);

			// Compute pairwise affinities
			elementPower(Y, 2, sqed);
			sumRows(sqed, sum_Y);

			// Transpose the sum vector.  The sumRows now
			// gratuitously tansposes sum_Y
			transpose(sum_Y);

			multAddTransB(-2.0, Y, Y, Ysqlmul);

			addRowVector(Ysqlmul, sum_Y);
			transpose(Ysqlmul);
			addRowVector(Ysqlmul, sum_Y);
			
			add(Ysqlmul, 1.0);
			divide(1.0,Ysqlmul);
			num.set(Ysqlmul);
			assignAtIndex(num, range(n), range(n), 0);
			divide(num , round(elementSum(num), DIGITS), Q);

			maximize(Q, 1e-12);
			// writeMatrix("Q", Q);
			
			// Compute gradient
			subtract(P, Q, L);
			elementMult(L, num);
			DMatrixRMaj rowsum = sumRows(L,null); // rowsum = nx1
			double [] rsum  = new double[rowsum.numRows];
			for (int i = 0; i < rsum.length; i++) {
				rsum[i] = rowsum.get(i,0);
			}
			setDiag(diag,rsum);
			subtract(diag, L, L);
			mult(L, Y, dY);
			scale(4.0, dY);
			
			// Perform the update
			if (iter < 20)
				momentum = initial_momentum;
			else
				momentum = final_momentum;
			
			boolean [][] boolMtrx = equal(biggerThan(dY,0.0),biggerThan(iY,0.0));
			
			
			setData(btNeg, abs(negate(boolMtrx)));
			setData(bt, abs(boolMtrx));
			
			DMatrixRMaj gainsSmall = new DMatrixRMaj(gains);
			DMatrixRMaj gainsBig   = new DMatrixRMaj(gains);
			add(gainsSmall,0.2);
			scale(0.8,gainsBig);
			
			elementMult(gainsSmall, btNeg);
			elementMult(gainsBig, bt);
			add(gainsSmall,gainsBig,gains);

			assignAllLessThan(gains, min_gain, min_gain);
			
			scale(momentum,iY);
			DMatrixRMaj gainsdY = new DMatrixRMaj(gains.numRows,dY.numCols);
			elementMult(gains , dY, gainsdY);
			scale(eta,gainsdY);
			subtractEquals(iY , gainsdY);
			addEquals(Y , iY);
			DMatrixRMaj colMeanY = colMean(Y, 0);
			DMatrixRMaj meanTile = tile(colMeanY, n, 1);
			subtractEquals(Y , meanTile);

			// Compute current value of cost function
			if (iter % 100 == 0)   {
				DMatrixRMaj Pdiv = new DMatrixRMaj(P);
				elementDiv(Pdiv , Q);
				elementLog(Pdiv,logdivide);
				replaceNaN(logdivide,Double.MIN_VALUE);
				elementMult(logdivide,P);
				replaceNaN(logdivide,Double.MIN_VALUE);
				double C = round(elementSum(logdivide), DIGITS);
				end = System.currentTimeMillis();
				 monitor.showMessage(TaskMonitor.Level.INFO,
				 				String.format("Iteration %d: error is %f (100 iterations in %4.2f seconds)\n", iter, C, (end - start) / 1000.0));
				if(C < 0) {
					monitor.showMessage(TaskMonitor.Level.WARN, "Warning: Error is negative, this is usually a very bad sign!");
				}
				start = System.currentTimeMillis();
			} /*else if(iter % 10 == 0) {
				end = System.currentTimeMillis();
				System.out.printf("Iteration %d: (10 iterations in %4.2f seconds)\n", iter, (end - start) / 1000.0);
				start = System.currentTimeMillis();
			} */

			// Stop lying about P-values
			if (iter == 100)
				divide(P , 4);
		}

		end = System.currentTimeMillis();
		monitor.showMessage(TaskMonitor.Level.INFO,
				String.format("Completed in %4.2f seconds)",(end-total)/1000.0));
		// Return solution
		return extractDoubleArray(Y);
	}
	
	public R Hbeta (double [][] D, double beta, int index){
		DMatrixRMaj P  = new DMatrixRMaj(D);
		if (index == 0)
			writeMatrix("P"+index, P);
		scale(-beta,P);
		if (index == 0)
			writeMatrix("HbetaP-scaled-"+index, P);
		elementExp(P,P);
		if (index == 0)
			writeMatrix("HbetaP-"+index, P);
		double sumP = elementSum(P);   // sumP confirmed scalar
			//System.out.println("sumP = "+sumP);
		DMatrixRMaj Dd  = new DMatrixRMaj(D);
		elementMult(Dd, P);
		double H = Math.log(sumP) + beta * elementSum(Dd) / sumP;
		scale(1/sumP,P);
		R r = new R();
		r.H = H;
		r.P = extractDoubleArray(P);
		return r;
	}

	public R x2p(double [][] X,double tol, double perplexity){
		int n               = X.length;
		double [][] sum_X   = sum(square(X), 1);
		writeMatrix("sum_X", sum_X);
		double [][] times   = scalarMult(times(X, mo.transpose(X)), -2);
		writeMatrix("times", times);
		double [][] prodSum = addColumnVector(mo.transpose(times), sum_X);
		writeMatrix("prodSum", prodSum);
		double [][] D       = MatrixOps.addRowVector(prodSum, mo.transpose(sum_X));
		// D seems correct at this point compared to Python version
		writeMatrix("D", D);
		double [][] P       = fillMatrix(n,n,0.0);
		double [] beta      = fillMatrix(n,n,1.0)[0];
		double logU         = Math.log(perplexity);
		// System.out.println("Starting x2p...");
		for (int i = 0; i < n; i++) {
			if (config.cancelled())
				break;
			if (i % 500 == 0)
				monitor.showMessage(TaskMonitor.Level.INFO, "Computing P-values for point " + i + " of " + n + "...");
			double betamin = Double.NEGATIVE_INFINITY;
			double betamax = Double.POSITIVE_INFINITY;
			double [][] Di = getValuesFromRow(D, i,concatenate(range(0,i),range(i+1,n)));
			if (i == 0)
				writeMatrix("Di", Di);

			R hbeta = Hbeta(Di, beta[i], i);
			double H = hbeta.H;
			double [][] thisP = hbeta.P;

			// Evaluate whether the perplexity is within tolerance
			double Hdiff = H - logU;
			int tries = 0;
			while(Math.abs(Hdiff) > tol && tries < 50){
				if (Hdiff > 0){
					betamin = beta[i];
					if (Double.isInfinite(betamax))
						beta[i] = beta[i] * 2;
					else 
						beta[i] = (beta[i] + betamax) / 2;
				} else{
					betamax = beta[i];
					if (Double.isInfinite(betamin))  
						beta[i] = beta[i] / 2;
					else 
						beta[i] = ( beta[i] + betamin) / 2;
				}

				hbeta = Hbeta(Di, beta[i], i);
				H = hbeta.H;
				thisP = hbeta.P;
				Hdiff = H - logU;
				tries = tries + 1;
			}
			assignValuesToRow(P, i,concatenate(range(0,i),range(i+1,n)),thisP[0]);
		}

		R r = new R();
		r.P = P;
		r.beta = beta;
		double sigma = mean(sqrt(scalarInverse(beta)));

		monitor.showMessage(TaskMonitor.Level.INFO, "Mean value of sigma: " + sigma);

		return r;
	}

	@Override
	public void abort() {
		abort = true;
	}

	public static void writeMatrix(String fileName, double[][] matrix) {
		String filePath = "/tmp/" + fileName;
		try{
			File file = new File(filePath);
			if(!file.exists()) {
				file.createNewFile();
			}
			PrintWriter writer = new PrintWriter(filePath, "UTF-8");
			writer.write(printMatrix(matrix));
			writer.close();
		}catch(IOException e){
			e.printStackTrace(System.out);
		}
	}
	
	public static String printMatrix(double[][] matrix) {
		StringBuilder sb = new StringBuilder();
		int nRows = matrix.length;
		int nColumns = matrix[0].length;
		sb.append("raw array("+nRows+", "+nColumns+")\n\t");
    for (int col = 0; col < nColumns; col++) {
      sb.append("null\t");
    }
    sb.append("\n");
    for (int row = 0; row < nRows; row++) {
      sb.append("null:\t"); //node.getIdentifier()
      for (int col = 0; col < nColumns; col++) {
        sb.append(""+matrix[row][col]+"\t");
      }
      sb.append("\n");
    }
    return sb.toString();

	}

	public static void writeMatrix(String fileName, double[] array) {
		String filePath = "/tmp/" + fileName;
		try{
			File file = new File(filePath);
			if(!file.exists()) {
				file.createNewFile();
			}
			PrintWriter writer = new PrintWriter(filePath, "UTF-8");
			writer.write(printMatrix(array));
			writer.close();
		}catch(IOException e){
			e.printStackTrace(System.out);
		}
	}
	
	public static String printMatrix(double[] array) {
		StringBuilder sb = new StringBuilder();
		int nRows = array.length;
		int nColumns = 1;
		sb.append("raw array("+nRows+", "+nColumns+")\n\t");
    sb.append("\n");
    for (int row = 0; row < nRows; row++) {
      sb.append("null:\t"+array[row]+"\n"); //node.getIdentifier()
    }
    return sb.toString();

	}
	
	public static void writeMatrix(String fileName, DMatrixRMaj matrix) {
		String filePath = "/tmp/" + fileName;
		try{
			File file = new File(filePath);
			if(!file.exists()) {
				file.createNewFile();
			}
			PrintWriter writer = new PrintWriter(filePath, "UTF-8");
			writer.write(printMatrix(matrix));
			writer.close();
		}catch(IOException e){
			e.printStackTrace(System.out);
		}
	}

	public static String printMatrix(DMatrixRMaj matrix) {
		StringBuilder sb = new StringBuilder();
		int nRows = matrix.getNumRows();
		int nColumns = matrix.getNumCols();
		sb.append("EJML Matrix("+nRows+", "+nColumns+")\n\t");
    for (int col = 0; col < nColumns; col++) {
      sb.append("null\t");
    }
    sb.append("\n");
    for (int row = 0; row < nRows; row++) {
      sb.append("null:\t"); //node.getIdentifier()
      for (int col = 0; col < nColumns; col++) {
        sb.append(""+matrix.get(row,col)+"\t");
      }
      sb.append("\n");
    }
    return sb.toString();

	}

	public double round(double value, int numberOfDigitsAfterDecimalPoint) {
		//System.out.println("value = "+value);
		BigDecimal bigDecimal = new BigDecimal(value);
		bigDecimal = bigDecimal.setScale(numberOfDigitsAfterDecimalPoint, BigDecimal.ROUND_HALF_UP);
		return bigDecimal.doubleValue();
	}


}
