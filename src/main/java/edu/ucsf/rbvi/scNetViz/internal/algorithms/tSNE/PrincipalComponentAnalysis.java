package edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE;
/*
 * Copyright (c) 2009-2014, Peter Abeles. All Rights Reserved.
 *
 * This file is part of Efficient Java Matrix Library (EJML).
 * 
 * Adapted by 2014, Leif Jonsson, added pca method.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */


/*
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.factory.DecompositionFactory_DDRM;
import org.ejml.interfaces.decomposition.SingularValueDecomposition_F64;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.NormOps_DDRM;
import org.ejml.dense.row.SingularOps_DDRM;
*/

import org.ojalgo.OjAlgoUtils;
import org.ojalgo.matrix.decomposition.DecompositionStore;
import org.ojalgo.matrix.decomposition.Eigenvalue;
import org.ojalgo.matrix.decomposition.SingularValue;
import org.ojalgo.function.aggregator.Aggregator;
import org.ojalgo.function.aggregator.AggregatorFunction;
import org.ojalgo.matrix.store.MatrixStore;
import org.ojalgo.matrix.store.PhysicalStore;
import org.ojalgo.matrix.store.PrimitiveDenseStore;
import org.ojalgo.matrix.store.RawStore;
import org.ojalgo.scalar.PrimitiveScalar;
import org.ojalgo.scalar.Scalar;


/**
 * <p>
 * The following is a simple example of how to perform basic principal component analysis in EJML.
 * </p>
 *
 * <p>
 * Principal Component Analysis (PCA) is typically used to develop a linear model for a set of data
 * (e.g. face images) which can then be used to test for membership.  PCA works by converting the
 * set of data to a new basis that is a subspace of the original set.  The subspace is selected
 * to maximize information.
 * </p>
 * <p>
 * PCA is typically derived as an eigenvalue problem.  However in this implementation {@link org.ejml.interfaces.decomposition.SingularValueDecomposition SVD}
 * is used instead because it will produce a more numerically stable solution.  Computation using EVD requires explicitly
 * computing the variance of each sample set. The variance is computed by squaring the residual, which can
 * cause loss of precision.
 * </p>
 *
 * <p>
 * Usage:<br>
 * 1) call setup()<br>
 * 2) For each sample (e.g. an image ) call addSample()<br>
 * 3) After all the samples have been added call computeBasis()<br>
 * 4) Call  sampleToEigenSpace() , eigenToSampleSpace() , errorMembership() , response()
 * </p>
 *
 * @author Peter Abeles
 */
public class PrincipalComponentAnalysis {

		protected final PhysicalStore.Factory<Double, PrimitiveDenseStore> storeFactory;
    // principal component subspace is stored in the rows
    private MatrixStore<Double> V_t;

    // how many principal components are used
    private int numComponents;

    // where the data is stored
    private PhysicalStore<Double> A;
    private int sampleIndex;

    // mean values of each element across all the samples
    double mean[];

    public PrincipalComponentAnalysis() {
			storeFactory = PrimitiveDenseStore.FACTORY;
    }

    /**
     * Must be called before any other functions. Declares and sets up internal data structures.
     *
     * @param numSamples Number of samples that will be processed.
     * @param sampleSize Number of elements in each sample.
     */
    public void setup( int numSamples , int sampleSize ) {
        mean = new double[ sampleSize ];
        A = storeFactory.makeZero(numSamples, sampleSize);
        sampleIndex = 0;
        numComponents = -1;
    }

    /**
     * Adds a new sample of the raw data to internal data structure for later processing.  All the samples
     * must be added before computeBasis is called.
     *
     * @param sampleData Sample from original raw data.
     */
    public void addSample( double[] sampleData ) {
        if( A.countColumns() != sampleData.length )
            throw new IllegalArgumentException("Unexpected sample size");
        if( sampleIndex >= A.countRows() )
            throw new IllegalArgumentException("Too many samples");

        for( int i = 0; i < sampleData.length; i++ ) {
            A.set(sampleIndex,i,sampleData[i]);
        }
        sampleIndex++;
    }

    /**
     * Computes a basis (the principal components) from the most dominant eigenvectors.
     *
     * @param numComponents Number of vectors it will use to describe the data.  Typically much
     * smaller than the number of elements in the input vector.
     */
    public void computeBasis( int numComponents ) {
        if( numComponents > A.countColumns() )
            throw new IllegalArgumentException("More components requested that the data's length.");
        if( sampleIndex != A.countRows() )
            throw new IllegalArgumentException("Not all the data has been added");
        if( numComponents > sampleIndex )
            throw new IllegalArgumentException("More data needed to compute the desired number of components");

        this.numComponents = numComponents;

				// System.out.println("Computing mean");

				// System.out.println("Centralizing");
        // Centralize the rows
        for( int row = 0; row < A.countRows(); row++ ) {
						final AggregatorFunction<Double> tmpVisitor = MySUM.get().reset();
						A.visitRow(row, 0L, tmpVisitor);
						mean[row] = tmpVisitor.doubleValue()/A.countColumns();

						for (int col = 0; col < A.countColumns(); col++) {
							double cell = A.get(row, col);
							A.set(row, col, cell - mean[row]);
            }
        }

				// System.out.println("Computing SVD");
				final SingularValue<Double> singularValue = SingularValue.make(A);
        singularValue.compute(A);
				// System.out.println("isOrdered = "+singularValue.isOrdered());

        V_t = singularValue.getQ2().transpose();
        // MatrixStore<Double> W = singularValue.getD();

				// System.out.println("sorting");

        // Singular values are in an arbitrary order initially
        // SingularOps_DDRM.descendingOrder(null,false,W,V_t,true);

				// System.out.println("reshaping");

        // strip off unneeded components and find the basis
        // V_t.reshape(numComponents,mean.length,true);
				// System.out.println("Limiting");
				V_t = V_t.logical().limits(numComponents, mean.length).get();
    }

    /**
     * Returns a vector from the PCA's basis.
     *
     * @param which Which component's vector is to be returned.
     * @return Vector from the PCA basis.
     *
    public double[] getBasisVector( int which ) {
        if( which < 0 || which >= numComponents )
            throw new IllegalArgumentException("Invalid component");

        DMatrixRMaj v = new DMatrixRMaj(1,A.numCols);
        CommonOps_DDRM.extract(V_t,which,which+1,0,A.numCols,v,0,0);

        return v.data;
    }
		*/

    /**
     * Converts a vector from sample space into eigen space.
     *
     * @param sampleData Sample space data.
     * @return Eigen space projection.
     */
    public double[] sampleToEigenSpace( double[] sampleData, int row ) {
        if( sampleData.length != A.countColumns() )
            throw new IllegalArgumentException("Unexpected sample length");

				MatrixStore<Double> mean = new RawStore(this.mean, 1);
        // DMatrixRMaj mean = DMatrixRMaj.wrap(A.getNumCols(),1,this.mean);

				MatrixStore<Double> s = new RawStore(sampleData, 1).copy();
        // DMatrixRMaj s = new DMatrixRMaj(A.getNumCols(),1,true,sampleData);

				// MatrixStore<Double> r = storeFactory.makeZero(numComponents, 1);
        // DMatrixRMaj r = new DMatrixRMaj(numComponents,1);

        // CommonOps_DDRM.subtract(s, mean, s);
				s = s.subtract(mean);

				// FastTSne.writeMatrix("s-"+row,s);

        // CommonOps_DDRM.mult(V_t,s,r);
				MatrixStore<Double> r = V_t.multiply(s);
				// FastTSne.writeMatrix("V_t-"+row,V_t);
				// FastTSne.writeMatrix("r-"+row,r);

        return r.toRawCopy1D();
    }

    /**
     * Converts a vector from eigen space into sample space.
     *
     * @param eigenData Eigen space data.
     * @return Sample space projection.
     *
    public double[] eigenToSampleSpace( double[] eigenData ) {
        if( eigenData.length != numComponents )
            throw new IllegalArgumentException("Unexpected sample length");

        DMatrixRMaj s = new DMatrixRMaj(A.getNumCols(),1);
        DMatrixRMaj r = DMatrixRMaj.wrap(numComponents,1,eigenData);
        
        CommonOps_DDRM.multTransA(V_t,r,s);

        DMatrixRMaj mean = DMatrixRMaj.wrap(A.getNumCols(),1,this.mean);
        CommonOps_DDRM.add(s,mean,s);

        return s.data;
    }
		*/


    /**
     * <p>
     * The membership error for a sample.  If the error is less than a threshold then
     * it can be considered a member.  The threshold's value depends on the data set.
     * </p>
     * <p>
     * The error is computed by projecting the sample into eigenspace then projecting
     * it back into sample space and
     * </p>
     * 
     * @param sampleA The sample whose membership status is being considered.
     * @return Its membership error.
     *
    public double errorMembership( double[] sampleA ) {
        double[] eig = sampleToEigenSpace(sampleA, 0);
        double[] reproj = eigenToSampleSpace(eig);


        double total = 0;
        for( int i = 0; i < reproj.length; i++ ) {
            double d = sampleA[i] - reproj[i];
            total += d*d;
        }

        return Math.sqrt(total);
    }
		*/

    /**
     * Computes the dot product of each basis vector against the sample.  Can be used as a measure
     * for membership in the training sample set.  High values correspond to a better fit.
     *
     * @param sample Sample of original data.
     * @return Higher value indicates it is more likely to be a member of input dataset.
     *
    public double response( double[] sample ) {
        if( sample.length != A.numCols )
            throw new IllegalArgumentException("Expected input vector to be in sample space");

        DMatrixRMaj dots = new DMatrixRMaj(numComponents,1);
        DMatrixRMaj s = DMatrixRMaj.wrap(A.numCols,1,sample);

        CommonOps_DDRM.mult(V_t,s,dots);

        return NormOps_DDRM.normF(dots);
    }
		*/
    
		// TODO: wrap all of this in a matrix rather than a 2-dimensional array
    public double [][] pca(double [][]matrix, int no_dims) {
			double [][] trafoed = new double[matrix.length][matrix[0].length];
			// System.out.println("setup");
			setup(matrix.length, matrix[0].length);
			// System.out.println("adding samples");
			for (int i = 0; i < matrix.length; i++) {
				addSample(matrix[i]);
			}
			// System.out.println("computing basis");
			computeBasis(no_dims);
			// System.out.println("Converting to eigenspace");
			for (int i = 0; i < matrix.length; i++) {
				trafoed[i] = sampleToEigenSpace(matrix[i], i);
				for (int j = 0; j < trafoed[i].length; j++) {
					trafoed[i][j] *= -1;
				}
			}
			// System.out.println("done");
			return trafoed;
    }

		public static final ThreadLocal<AggregatorFunction<Double>> MySUM = new ThreadLocal<AggregatorFunction<Double>>() {

    @Override
    protected AggregatorFunction<Double> initialValue() {
      return new AggregatorFunction<Double>() {

        private double sum = 0.0;

        public void invoke(final Double anArg) {
          if (anArg != null) invoke(anArg.doubleValue());
        }
        public void invoke(final double anArg) {
          if (!Double.isNaN(anArg))
            sum += anArg;
        }

        public double doubleValue() { return sum; }
        public Scalar<Double> toScalar() {
          return PrimitiveScalar.of(this.doubleValue());
        }
        public AggregatorFunction<Double> reset() { sum = 0.0; return this; }
        public void merge(final Double result) {
          this.invoke(result.doubleValue());
        }
        public Double merge(final Double result1, final Double result2) {
          return result1 + result2;
        }
        /*
        public Double getNumber() {
          return Double.valueOf(this.doubleValue());
        }
        */
        public Double get() {
          return Double.valueOf(this.doubleValue());
        }
      };
    }
  };

}
