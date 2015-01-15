#include <EventStruct.h>
#include <MovingParticle.h>
#include <mex.h>

EventType int2EventType(int i)
{
	EventType type = UnknownEvent;
	switch (i)
	{
	case 0:
		break;
	case 1:
		type = CollisionEvent;
		break;
	case 2:
		type = SplitEvent;
		break;
	case 3:
		type = MergeEvent;
		break;
	default:
		break;
	}
	return type;
}


EventStruct::EventStruct(float t, EventType type, const MovingParticle* p, const MovingParticle* q)
{
	this->t = t;
	this->type = type;
	this->p = p;
	this->q = q;
	if (q == NULL)
	{
		r = NULL;
	}
	else
	{
		r = q->getNext();
	}
}

void EventStruct::print()
{
	CParticleF pp = p->getP();
	if (type == SplitEvent)
	{
		CParticleF qp = q->getP();
		CParticleF rp = r->getP();
		MovingParticle* c1 = p->getChildren(0);
		MovingParticle* c2 = p->getChildren(1);
		printf("event>> Split @ %f: %d(%3.3f,%3.3f)->%d(%3.3f,%3.3f), %d(%3.3f,%3.3f) => %d, %d\n",
			t, p->getId(), pp.m_X, pp.m_Y, q->getId(), qp.m_X, qp.m_Y, r->getId(), rp.m_X, rp.m_Y,
			c1==NULL ? -1: c1->getId(), c2==NULL ? -1: c2->getId());
	}
	else if (type == CollisionEvent)
	{
		CParticleF qp = q->getP();
		MovingParticle* c1 = p->getChildren(0);
		printf("event>> Collision @ %f: %d(%3.3f,%3.3f)->%d(%3.3f,%3.3f) => %d\n",
			t, p->getId(), pp.m_X, pp.m_Y, q->getId(), qp.m_X, qp.m_Y,
			c1==NULL ? -1: c1->getId());
	}
	else
	{
		printf("event>> Unknown @ %f: %d(%3.3f,%3.3f)\n",
			t, p->getId(), pp.m_X, pp.m_Y);
	}
}
