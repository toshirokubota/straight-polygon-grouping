#pragma once
#include <limits>
class MovingParticle;

/*
CollisionEvent: two adjacent vertices collide and replaced by a new one.
SplitEvent: a concave vertex collide with a side and split the polygon into two.
MergeEvent: a vertex collide with a side that is parallel to the colliding vertex.
*/
enum EventType { UnknownEvent, CollisionEvent, SplitEvent, MergeEvent };

EventType int2EventType(int i);

struct EventStruct
{
	EventStruct(float t = std::numeric_limits<float>::infinity(), EventType type = UnknownEvent,	
				const MovingParticle* p = NULL, 
				const MovingParticle* q = NULL);

	bool operator <(EventStruct& ev)
	{
		return t < ev.t;
	}
	float t;
	EventType type;
	const MovingParticle* p;
	const MovingParticle* q;
	//r is different for each event type.
	//It is NULL for a Collision event.
	//It is q->next for a Split event.
	//It is p->prev for a Merge event if p->prev is in the same merge cluster. It is NULL, otherwise.
	const MovingParticle* r;

	void print();
};


