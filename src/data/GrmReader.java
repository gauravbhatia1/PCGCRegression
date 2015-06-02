package data;

import java.io.BufferedInputStream;
import java.io.EOFException;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;

import com.google.common.io.LittleEndianDataInputStream;



public class GrmReader {
	protected int m_index1, m_index2;
	private int m_next_index1, m_next_index2;
	private float m_value;
	private float m_next_value;
	private boolean m_has_next = true;
	private String m_bin_file_name;
	private String m_id_file_name;
	private LittleEndianDataInputStream m_input_stream;
	private IndividualVariables m_ids;
	
	public GrmReader(String bin_filename, String id_filename) throws FileNotFoundException
	{
		m_index1 = -1;
		m_index2 = -1;
		m_next_index1 = -1;
		m_next_index2 = -1;
		m_value = -1;
		m_bin_file_name = bin_filename;
		m_id_file_name = id_filename;
		m_input_stream = new LittleEndianDataInputStream(new BufferedInputStream(new FileInputStream(m_bin_file_name)));
		m_ids = IndividualVariables.read(m_id_file_name);
		next();
	}
	
	public int getIndex1()
	{
		return m_index1;
	}
	
	public int getIndex2()
	{
		return m_index2;
	}
	
	public float getValue()
	{
		return m_value;
	}
	
	public void next()
	{
		if( !m_has_next )
		{
			throw new RuntimeException();
		}
		
		m_value = m_next_value;
		m_index1 = m_next_index1;
		m_index2 = m_next_index2;
		try {
			m_next_value = m_input_stream.readFloat();
			m_next_index2++;
			if( m_next_index2 > m_next_index1 )
			{
				m_next_index1++;
				m_next_index2 = 0;
			}
		} catch (EOFException e) {
			m_has_next = false;
		}catch (IOException e) {
			throw new RuntimeException(e);
		}
	}
	
	public boolean hasNext() throws IOException
	{
		return m_has_next;
	}

	public IndividualVariables getIds() {
		return m_ids;
	}
	
	
}
