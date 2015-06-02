package data;

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;

import com.google.common.io.LittleEndianDataInputStream;


public class Grm {
	private double[][] m_grm;
	private Map<IndividualId, Integer> m_individual_id_map;
	
	private Grm(double[][] grm, Map<IndividualId, Integer> id_map)
	{
		m_grm = grm;
		m_individual_id_map = id_map;
	}
	
	public double getCovariance(IndividualId id1, IndividualId id2)
	{
		Integer index1 = m_individual_id_map.get(id1);
		if( index1 == null )
		{
			throw new RuntimeException("no individual with id " + id1);
		}
		Integer index2 = m_individual_id_map.get(id2);
		if( index2 == null )
		{
			throw new RuntimeException("no individual with id " + id2);
		}
		return m_grm[index1][index2];
	}
	
	public double getCovariance(int index1, int index2)
	{
		if( index1 < 0 || index1 >= m_grm.length )
		{
			throw new RuntimeException(index1 + " is out of range");
		}
		if( index2 < 0 || index2 >= m_grm.length )
		{
			throw new RuntimeException(index2 + " is out of range");
		}
		return m_grm[index1][index2];
	}
	
	public double getQuickCovariance(int index1, int index2) {
		return m_grm[index1][index2];
	}
	
	public static Grm readGrmBin(String prefix) throws IOException
	{
		Map<IndividualId, Integer> indexes = new HashMap<IndividualId, Integer>();
		String id_file = prefix + ".grm.id";
		Scanner scanner = new Scanner(new BufferedInputStream(new FileInputStream(id_file)));
		int index = 0;
		while(scanner.hasNext())
		{
			String fid = scanner.next();
			String iid = scanner.next();
			IndividualId id = new IndividualId(fid, iid);
			indexes.put(id, index);
			index++;
		}
		scanner.close();
		int indvidual_count = index;
		double[][] grm = new double[indvidual_count][];
		for(int i = 0; i < indvidual_count; i++)
		{
			grm[i] = new double[indvidual_count];
		}
		String binfile = prefix + ".grm.bin";
		LittleEndianDataInputStream dis = new LittleEndianDataInputStream(new BufferedInputStream(new FileInputStream(binfile)));
		int i = 0;
		int j = 0;
		while(dis.available() > 0)
		{
			float value = dis.readFloat();
			grm[i][j] = value;
			grm[j][i] = value;
			j++;
			if( j > i )
			{
				j = 0;
				i++;
			}
		}
		dis.close();
		return new Grm(grm, indexes);
	}

	public Set<IndividualId> getIndividualIds() {
		return m_individual_id_map.keySet();
	}
	
}
