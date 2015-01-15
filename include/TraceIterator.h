#ifndef ___TRACE_ITERATOR_H___
#define ___TRACE_ITERATOR_H___

struct TraceIterator
{
	TraceIterator(int x, int y, int px, int py)
	{
		int dx = x - px;
		int dy = y - py;
		int ix = dx > 0 ? 1: (dx < 0 ? -1: 0);
		int iy = dy > 0 ? 1: (dy < 0 ? -1: 0);
		if(ix == 0 && iy==0) index = 0;
		else
		{
			if(ix==-1 && iy==-1) index = 0;
			else if(ix==0 && iy==-1) index = 1;
			else if(ix==1 && iy==-1) index = 2;
			else if(ix==-1 && iy==0) index = 7;
			else if(ix==1 && iy==0) index = 3;
			else if(ix==-1 && iy==1) index = 6;
			else if(ix==0 && iy==1) index = 5;
			else if(ix==1 && iy==1) index = 4;
		}
		_XOffset8[0]=_XOffset8[6]=_XOffset8[7]=-1;
		_XOffset8[1]=_XOffset8[5]=0;
		_XOffset8[2]=_XOffset8[3]=_XOffset8[4]=1;
		_YOffset8[0]=_YOffset8[1]=_YOffset8[2]=-1;
		_YOffset8[3]=_YOffset8[7]=0;
		_YOffset8[4]=_YOffset8[5]=_YOffset8[6]=1;
		count = 0;
	}
	TraceIterator(int x, int y)
	{
		*this = TraceIterator(x, y, x, y);
	}
	TraceIterator(const CParticleF& q, const CParticleF& p)
	{
		*this = TraceIterator(q.m_X, q.m_Y, p.m_X, p.m_Y);
	}
	bool Next(int& x, int& y)
	{
		if(count >= 8) return false; 
		else
		{
			int k;
			if(count % 2 == 1)
			{
				k = (index - (count+1)/2 + 8) % 8;
			}
			else
			{
				k = (index + count/2) % 8;
			}
			x = _XOffset8[k];
			y = _YOffset8[k];
			count ++;
			return true;
		}
	}
	int XOffset(int k) {return _XOffset8[k];}
	int YOffset(int k) {return _YOffset8[k];}
private:
	int index;
	int count;
	int _XOffset8[8];
	int _YOffset8[8];
};

#endif /* ___TRACE_ITERATOR_H___ */