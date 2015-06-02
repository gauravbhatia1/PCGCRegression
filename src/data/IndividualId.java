package data;

public class IndividualId {
	private String m_fid, m_iid;

	public IndividualId(String fid, String iid) {
		super();
		m_fid = fid;
		m_iid = iid;
	}

	public String getFid() {
		return m_fid;
	}

	public String getIid() {
		return m_iid;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((m_fid == null) ? 0 : m_fid.hashCode());
		result = prime * result + ((m_iid == null) ? 0 : m_iid.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		IndividualId other = (IndividualId) obj;
		if (m_fid == null) {
			if (other.m_fid != null)
				return false;
		} else if (!m_fid.equals(other.m_fid))
			return false;
		if (m_iid == null) {
			if (other.m_iid != null)
				return false;
		} else if (!m_iid.equals(other.m_iid))
			return false;
		return true;
	}

	@Override
	public String toString() {
		return "IndividualId [m_fid=" + m_fid + ", m_iid=" + m_iid + "]";
	}
}
