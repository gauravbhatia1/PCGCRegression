package estimation;

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.Scanner;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.jet.math.Functions;
import data.GrmReader;
import data.IndividualId;
import data.IndividualVariables;

public class LowMemMultipleHeRegressionEstimatorBin {
	private static final String GRMS = "grmlist";
	private static final String PHENOS = "phenos";
	private static final String MULTIPLIER = "multiplier";
	private static final String TOTALVAR = "totalvar";
	private static final String OUTFILE = "outfile";
	
	public static void main(String[] args) throws IOException
	{
		String paramfile = args[0];
		Properties p = new Properties();
		p.load(new BufferedInputStream(new FileInputStream(paramfile)));
		String grmlist = p.getProperty(GRMS);
		String phenofile = p.getProperty(PHENOS);
		String multiplier_file = p.getProperty(MULTIPLIER);
		double total_var = Double.parseDouble(p.getProperty(TOTALVAR));
		
		List<GrmReader> l = getGrms(grmlist);
		IndividualVariables ids = l.get(0).getIds();
		List<IndividualId> ordered_ids = ids.getOrderedIndividualIds();
		IndividualVariables pheno_set = IndividualVariables.read(phenofile);
		IndividualVariables multiplier_set = IndividualVariables.read(multiplier_file);
		List<IndividualId> filtered_ids = new ArrayList<IndividualId>();
		Map<Integer, Integer> mapped = produceFinal(ordered_ids, pheno_set, multiplier_set);
		for(int i = 0; i< ordered_ids.size(); i++)
		{
			if( (mapped == null) || (mapped.containsKey(i)))
			{
				filtered_ids.add(ordered_ids.get(i));
			}
		}
		
		double[] pheno_array = pheno_set.getArray(filtered_ids, 0);
		double[] mult_array = multiplier_set.getArray(filtered_ids, 0);
		DoubleMatrix2D den = new DenseDoubleMatrix2D(l.size(), l.size() );
		DoubleMatrix2D num = new DenseDoubleMatrix2D(l.size(), 1);
		int count = 0;
		while(l.get(0).hasNext())
		{
			count += addNext(l, den, num, pheno_array, mult_array, mapped);
		}
		if( p.getProperty(OUTFILE) != null )
		{
			String outfile = p.getProperty(OUTFILE);
			ObjectOutputStream oos = new ObjectOutputStream(new FileOutputStream(outfile));
			oos.writeInt(count);
			oos.writeObject(den);
			oos.writeObject(num);
			oos.close();
		}
		Algebra a = new Algebra();
		DoubleMatrix2D inv = a.inverse(den);
		DoubleMatrix2D coeffs = a.mult(inv, num);
		for(int i = 0; i < coeffs.rows(); i++)
		{
			double vg = coeffs.getQuick(i, 0);
			vg /= (1 + total_var);
			System.out.println("V(G" + i + ")/Vp:\t" + vg);
		}
	}
	
	private static Map<Integer, Integer> produceFinal(List<IndividualId> ordered_grm_id, IndividualVariables phenotypes, IndividualVariables covariates)
	{
		Map<Integer, Integer> mapped = new HashMap<Integer, Integer>();
		int index = -1;
		int mapped_index = -1;
		for(IndividualId id : ordered_grm_id)
		{
			index++;
			if( !phenotypes.contains(id) )
			{
				continue;
			}
			if( !covariates.contains(id) )
			{
				continue;
			}
			mapped_index++;
			mapped.put(index, mapped_index);
		}
		//all ids were kept, return null
		if( mapped.size() == ordered_grm_id.size() )
		{
			return null;
		}
		return mapped;
	}

	private static int addNext(List<GrmReader> l, DoubleMatrix2D den, DoubleMatrix2D num, double[] phenos, double[] multiplier, Map<Integer,Integer> mapped_indexes) throws IOException {
		double[][] vals = new double[l.size()][];
		int index1 = 0, index2 = 0;
		boolean keep_looking = true;
		while(keep_looking)
		{
			if(!l.get(0).hasNext())
			{
				return 0;
			}
			l.get(0).next();
			vals[0] = new double[1];
			vals[0][0] = l.get(0).getValue();
			index1 = l.get(0).getIndex1();
			index2 = l.get(0).getIndex2();
			for(int i = 1; i < l.size(); i++)
			{
				GrmReader s  = l.get(i);
				if(!s.hasNext())
				{
					throw new RuntimeException();
				}
				s.next();
				vals[i] = new double[1];
				vals[i][0] = s.getValue();
			}
			keep_looking = false;
			if(index1 == index2)
			{
				keep_looking = true;
			}
			else if(mapped_indexes != null) 
			{
				if( !mapped_indexes.containsKey(index1) )
				{
					keep_looking = true;
				}
				else if( !mapped_indexes.containsKey(index2) )
				{
					keep_looking = true;
				}
			}
		}
		int mapped_index1 = index1;
		int mapped_index2 = index2;
		if( mapped_indexes != null )
		{
			mapped_index1 = mapped_indexes.get(index1);
			mapped_index2 = mapped_indexes.get(index2);
		}
		if( multiplier != null )
		{
			for(int index = 0; index < vals.length; index++)
			{
				vals[index][0] *= multiplier[mapped_index1]*multiplier[mapped_index2];
			}
		}
		double yval = phenos[mapped_index1]*phenos[mapped_index2];
//		if( s_count < 10 )
//		{
//			StringBuffer buf = new StringBuffer();
//			for(int index = 1; index < vals.length; index++)
//			{
//				buf.append(vals[index][0]);
//				buf.append('\t');
//			}
//			buf.append(yval);
//			System.out.println(buf);
//			s_count++;
//		}
		DoubleMatrix2D d2d = new DenseDoubleMatrix2D(vals);
		d2d.zMult(d2d, den, 1, 1, false, true);
		num.assign(d2d, Functions.plusMult(yval));
		return 1;
	}


	public static List<GrmReader> getGrms(String grmlist) throws IOException
	{
		Scanner scanner = null;
		IndividualVariables ids = null;
		List<GrmReader> l;
		try {
			scanner = new Scanner(new BufferedInputStream(new FileInputStream(grmlist)));
			l = new ArrayList<GrmReader>();
			while(scanner.hasNext())
			{
				String prefix = scanner.next();
				String idfile = prefix + ".grm.id";
				String binfile = prefix + ".grm.bin";
				GrmReader reader = new GrmReader(binfile, idfile);
				if( ids == null )
				{
					ids = reader.getIds();
				}
				else if( !reader.getIds().getOrderedIndividualIds().equals(ids.getOrderedIndividualIds()) )
				{
					throw new RuntimeException(idfile + " not identical to previous\n");
				}
				l.add(reader);
			}
		} finally
		{
			if( scanner != null )
			{
				scanner.close();
			}
		}
		return l;
	}
	
}
