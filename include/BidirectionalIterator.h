#ifndef ___BIDIRECTIONAL_ITERATOR_H___
#define ___BIDIRECTIONAL_ITERATOR_H___
#include <vector>

struct BidirectionalIterator
{
public:
	BidirectionalIterator(int source, int size, int target=-1)
	{
		_source = source;
		_size = size;
		_forward = true;
		_current = 0;
		_end = false;
		if(target<0) _target = _source;
		else _target = target;
		for(int i=source; i<size; ++i)
		{
			_order.push_back(i);
		}
		for(int i=0; i<source; i++)
		{
			_order.push_back(i);
		}
		for(int i=target; i>=0; i--)
		{
			_reverse.push_back(i);
		}
		for(int i=size-1; i>target; i--)
		{
			_reverse.push_back(i);
		}

		/*for(int i=source; i<size; ++i)
		{
			_order.push_back(i);
		}
		for(int i=size-1; i>=0; i--)
		{
			_order.push_back(i);
		}
		for(int i=0; i<source; ++i)
		{
			_order.push_back(i);
		}*/
		
		//regular order
		/*for(int i=0; i<size; ++i)
		{
			_order.push_back(i);
			_reverse.push_back(size-1-i);
		}*/

	}
	int Restart() 
	{
		_current = 0;
		_end = false;
		if(_forward)
			return _order[_current];
		else
			return _reverse[_current];
	}
	int Next()
	{
		if(_forward) 
		{
			if(_current < _order.size() - 1)
			{
				_current++;
			}
			else
			{
				_end = true;
			} 
			return _order[_current];
		}
		else
		{
			if(_current < _reverse.size()-1)
			{
				_current++;
			}
			else
			{
				_end = true;
			}
			return _reverse[_current];
		}
	}
	bool End()
	{
		return _end;
	}
	void Turn()
	{
		_forward = _forward ? false : true;
	}
	bool Active(float afrom, float ato)
	{
		if(_forward)
		{
			return afrom <= ato;
		}
		else
		{
			return afrom > ato;
		}
	}

private:
	int _source;
	int _target;
	int _size;
	int _current;
	bool _forward;
	bool _end;
	std::vector<int> _order;
	std::vector<int> _reverse;
};

#endif /* ___BIDIRECTIONAL_ITERATOR_H___ */