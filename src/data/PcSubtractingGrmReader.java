package data;

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

public class PcSubtractingGrmReader extends GrmReader {
	private double[][] m_scaled_vecs;
	
	public PcSubtractingGrmReader(String bin_filename, String id_filename, String vec_file, String val_file) throws FileNotFoundException {
		super(bin_filename, id_filename);
		IndividualVariables vecs = IndividualVariables.read(vec_file);
		double[] vals = readVals(val_file);
		List<IndividualId> ids = getIds().getOrderedIndividualIds();
		int variable_count = vecs.getVariableCount();
		m_scaled_vecs = new double[variable_count][];
		for(int vec_index = 0; vec_index < variable_count; vec_index++)
		{
			double[] array = vecs.getArray(ids, vec_index);
			double scale = Math.sqrt(vals[vec_index]);
			for(int individual_index = 0; individual_index < array.length; individual_index++)
			{
				array[individual_index] *= scale;
			}
			m_scaled_vecs[vec_index] = array;
		}
	}

	private static double[] readVals(String val_file) throws FileNotFoundException {
		List<Double> vals = new ArrayList<Double>();
		Scanner scanner = new Scanner(new BufferedInputStream(new FileInputStream(val_file)));
		while(scanner.hasNext())
		{
			vals.add(scanner.nextDouble());
		}
		double[] val_array = new double[vals.size()];
		for(int i = 0; i < vals.size(); i++)
		{
			val_array[i] = vals.get(i);
		}
		return val_array;
	}

	@Override
	public float getValue() {
		float unsubtracted = super.getValue();
		float subtraction = 0;
		for(int vec_index = 0; vec_index < m_scaled_vecs.length; vec_index++)
		{
			subtraction += m_scaled_vecs[vec_index][m_index1]*m_scaled_vecs[vec_index][m_index2];
		}
		return (unsubtracted - subtraction);
	}
	
	

}
