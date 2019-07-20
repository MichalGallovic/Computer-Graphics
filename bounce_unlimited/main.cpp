#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <algorithm>
#include <vector>
#include <set>
#include <map>
#include <float.h>
#include <math.h>
#include <assert.h>
#include <cstdlib>
#include <cstring>
#include <string>
#include <sstream>

using namespace std;
/**
 SINGLE PRECISION
*/
#ifndef CYCLONE_PRECISION_H
#define CYCLONE_PRECISION_H



namespace cyclone {

typedef float real;
    
    /** Defines the highest value for the real number. */
#define REAL_MAX FLT_MAX
    
    /** Defines the precision of the square root operator. */
#define real_sqrt sqrtf
    /** Defines the precision of the absolute magnitude operator. */
#define real_abs fabsf
    /** Defines the precision of the sine operator. */
#define real_sin sinf
    
    /** Defines the precision of the cosine operator. */
#define real_cos cosf
    
    /** Defines the precision of the exponent operator. */
#define real_exp expf
    /** Defines the precision of the power operator. */
#define real_pow powf
    
    /** Defines the precision of the floating point modulo operator. */
#define real_fmod fmodf
    
#define R_PI 3.14159f

}

#endif // PRECISION_H

// VEKTOROVE OPERACIE
#ifndef CYCLONE_CORE_H
#define CYCLONE_CORE_H

namespace cyclone {
    
    class Vector3
    {
    public:
        /** Holds the value along the x axis. */
        real x;
        
        /** Holds the value along the y axis. */
        real y;
        
        /** Holds the value along the z axis. */
        real z;
        
    private:
        /** Padding to ensure 4 word alignment. */
        real pad;
        
    public:
        /** The default constructor creates a zero vector. */
        Vector3() : x(0), y(0), z(0) {}
        
        /**
         * The explicit constructor creates a vector with the given
         * components.
         */
        Vector3(const real x, const real y, const real z)
        : x(x), y(y), z(z) {}
        
        const static Vector3 GRAVITY;
        const static Vector3 HIGH_GRAVITY;
        const static Vector3 UP;
        const static Vector3 RIGHT;
        const static Vector3 OUT_OF_SCREEN;
        const static Vector3 X;
        const static Vector3 Y;
        const static Vector3 Z;
        
        // ... Other Vector3 code as before ...
        
        
        real operator[](unsigned i) const
        {
            if (i == 0) return x;
            if (i == 1) return y;
            return z;
        }
        
        real& operator[](unsigned i)
        {
            if (i == 0) return x;
            if (i == 1) return y;
            return z;
        }
        
        /** Adds the given vector to this. */
        void operator+=(const Vector3& v)
        {
            x += v.x;
            y += v.y;
            z += v.z;
        }
        
        /**
         * Returns the value of the given vector added to this.
         */
        Vector3 operator+(const Vector3& v) const
        {
            return Vector3(x+v.x, y+v.y, z+v.z);
        }
        
        /** Subtracts the given vector from this. */
        void operator-=(const Vector3& v)
        {
            x -= v.x;
            y -= v.y;
            z -= v.z;
        }
        
        /**
         * Returns the value of the given vector subtracted from this.
         */
        Vector3 operator-(const Vector3& v) const
        {
            return Vector3(x-v.x, y-v.y, z-v.z);
        }
        
        /** Multiplies this vector by the given scalar. */
        void operator*=(const real value)
        {
            x *= value;
            y *= value;
            z *= value;
        }
        
        /** Returns a copy of this vector scaled the given value. */
        Vector3 operator*(const real value) const
        {
            return Vector3(x*value, y*value, z*value);
        }
        
        /**
         * Calculates and returns a component-wise product of this
         * vector with the given vector.
         */
        Vector3 componentProduct(const Vector3 &vector) const
        {
            return Vector3(x * vector.x, y * vector.y, z * vector.z);
        }
        
        /**
         * Performs a component-wise product with the given vector and
         * sets this vector to its result.
         */
        void componentProductUpdate(const Vector3 &vector)
        {
            x *= vector.x;
            y *= vector.y;
            z *= vector.z;
        }
        
        /**
         * Calculates and returns the vector product of this vector
         * with the given vector.
         */
        Vector3 vectorProduct(const Vector3 &vector) const
        {
            return Vector3(y*vector.z-z*vector.y,
                           z*vector.x-x*vector.z,
                           x*vector.y-y*vector.x);
        }
        
        /**
         * Updates this vector to be the vector product of its current
         * value and the given vector.
         */
        void operator %=(const Vector3 &vector)
        {
            *this = vectorProduct(vector);
        }
        
        /**
         * Calculates and returns the vector product of this vector
         * with the given vector.
         */
        Vector3 operator%(const Vector3 &vector) const
        {
            return Vector3(y*vector.z-z*vector.y,
                           z*vector.x-x*vector.z,
                           x*vector.y-y*vector.x);
        }
        
        /**
         * Calculates and returns the scalar product of this vector
         * with the given vector.
         */
        real scalarProduct(const Vector3 &vector) const
        {
            return x*vector.x + y*vector.y + z*vector.z;
        }
        
        /**
         * Calculates and returns the scalar product of this vector
         * with the given vector.
         */
        real operator *(const Vector3 &vector) const
        {
            return x*vector.x + y*vector.y + z*vector.z;
        }
        
        /**
         * Adds the given vector to this, scaled by the given amount.
         */
        void addScaledVector(const Vector3& vector, real scale)
        {
            x += vector.x * scale;
            y += vector.y * scale;
            z += vector.z * scale;
        }
        
        /** Gets the magnitude of this vector. */
        real magnitude() const
        {
            return real_sqrt(x*x+y*y+z*z);
        }
        
        /** Gets the squared magnitude of this vector. */
        real squareMagnitude() const
        {
            return x*x+y*y+z*z;
        }
        
        /** Limits the size of the vector to the given maximum. */
        void trim(real size)
        {
            if (squareMagnitude() > size*size)
            {
                normalise();
                x *= size;
                y *= size;
                z *= size;
            }
        }
        
        /** Turns a non-zero vector into a vector of unit length. */
        void normalise()
        {
            real l = magnitude();
            if (l > 0)
            {
                (*this) *= ((real)1)/l;
            }
        }
        
        /** Returns the normalised version of a vector. */
        Vector3 unit() const
        {
            Vector3 result = *this;
            result.normalise();
            return result;
        }
        
        /** Checks if the two vectors have identical components. */
        bool operator==(const Vector3& other) const
        {
            return x == other.x &&
            y == other.y &&
            z == other.z;
        }
        
        /** Checks if the two vectors have non-identical components. */
        bool operator!=(const Vector3& other) const
        {
            return !(*this == other);
        }
        
        /**
         * Checks if this vector is component-by-component less than
         * the other.
         *
         * @note This does not behave like a single-value comparison:
         * !(a < b) does not imply (b >= a).
         */
        bool operator<(const Vector3& other) const
        {
            return x < other.x && y < other.y && z < other.z;
        }
        
        /**
         * Checks if this vector is component-by-component less than
         * the other.
         *
         * @note This does not behave like a single-value comparison:
         * !(a < b) does not imply (b >= a).
         */
        bool operator>(const Vector3& other) const
        {
            return x > other.x && y > other.y && z > other.z;
        }
        
        /**
         * Checks if this vector is component-by-component less than
         * the other.
         *
         * @note This does not behave like a single-value comparison:
         * !(a <= b) does not imply (b > a).
         */
        bool operator<=(const Vector3& other) const
        {
            return x <= other.x && y <= other.y && z <= other.z;
        }
        
        /**
         * Checks if this vector is component-by-component less than
         * the other.
         *
         * @note This does not behave like a single-value comparison:
         * !(a <= b) does not imply (b > a).
         */
        bool operator>=(const Vector3& other) const
        {
            return x >= other.x && y >= other.y && z >= other.z;
        }
        
        /** Zero all the components of the vector. */
        void clear()
        {
            x = y = z = 0;
        }
        
        /** Flips all the components of the vector. */
        void invert()
        {
            x = -x;
            y = -y;
            z = -z;
        }
        
    };
    
}

#endif // CORE_H


// CASTICE
#ifndef CYCLONE_PARTICLE_H
#define CYCLONE_PARTICLE_H

namespace cyclone {
    

    class Particle
    {
    public:
        
        
    protected:
       
        real inverseMass;

        real damping;
        
        Vector3 position;

        Vector3 velocity;
        
        Vector3 forceAccum;
        
        Vector3 acceleration;
        
        real size = 0.1f;
    public:
  
        void integrate(real duration);
        
 
        void setMass(const real mass);
  
        real getMass() const;
 
        void setInverseMass(const real inverseMass);
        
        real getInverseMass() const;
 
        bool hasFiniteMass() const;
        
        void setDamping(const real damping);
  
        real getDamping() const;
        
        void setPosition(const Vector3 &position);
        
        void setPosition(const real x, const real y, const real z);
        
        void getPosition(Vector3 *position) const;
        
        Vector3 getPosition() const;

        void setVelocity(const Vector3 &velocity);
        
        void setVelocity(const real x, const real y, const real z);

        void getVelocity(Vector3 *velocity) const;

        Vector3 getVelocity() const;

        void setAcceleration(const Vector3 &acceleration);
        
        void setAcceleration(const real x, const real y, const real z);
        
        void getAcceleration(Vector3 *acceleration) const;

        Vector3 getAcceleration() const;
        
        void clearAccumulator();

        void addForce(const Vector3 &force);
        
        void setSize(const real size) {
            this->size = size;
        }
        
        real getSize() {
            return size;
        }

    };
}

#endif // BODY_H


//GENERATORY SIL
#ifndef CYCLONE_PFGEN_H
#define CYCLONE_PFGEN_H

namespace cyclone {
    

    class ParticleForceGenerator
    {
    public:

        virtual void updateForce(Particle *particle, real duration) = 0;
    };
    
    class ParticleAnchoredSpring : public ParticleForceGenerator
    {
    protected:
        Vector3 *anchor;
        
        real springConstant;
    
        real restLength;
        
    public:
        ParticleAnchoredSpring();
        
        ParticleAnchoredSpring(Vector3 *anchor,
                               real springConstant,
                               real restLength);
        

        const Vector3* getAnchor() const { return anchor; }
        

        void init(Vector3 *anchor,
                  real springConstant,
                  real restLength);
        

        virtual void updateForce(Particle *particle, real duration);
    };
    
    class ParticleForceRegistry
    {
    protected:
        

        struct ParticleForceRegistration
        {
            Particle *particle;
            ParticleForceGenerator *fg;
        };

        typedef std::vector<ParticleForceRegistration> Registry;
        Registry registrations;
        
    public:

        void add(Particle* particle, ParticleForceGenerator *fg);

        void remove(Particle* particle, ParticleForceGenerator *fg);

        void clear();

        void updateForces(real duration);
    };
}

#endif // PARTICLE_FORCE_GENERATORS_H


// KONTAKTY MEDZI CASTICAMI
#ifndef CYCLONE_PCONTACTS_H
#define CYCLONE_PCONTACTS_H

namespace cyclone {
    
    class ParticleContactResolver;

    class ParticleContact
    {

        friend class ParticleContactResolver;
        
        
    public:

        Particle* particle[2];

        real restitution;

        Vector3 contactNormal;

        real penetration;
        
        Vector3 particleMovement[2];
        
    protected:

        void resolve(real duration);

        real calculateSeparatingVelocity() const;
        
    private:

        void resolveVelocity(real duration);
        
        void resolveInterpenetration(real duration);
        
    };
    
    class ParticleContactResolver
    {
    protected:

        unsigned iterations;

        unsigned iterationsUsed;
        
    public:

        ParticleContactResolver(unsigned iterations);
        
        void setIterations(unsigned iterations);
        
        void resolveContacts(ParticleContact *contactArray,
                             unsigned numContacts,
                             real duration);
    };
    
    class ParticleContactGenerator
    {
    public:

        virtual unsigned addContact(ParticleContact *contact,
                                    unsigned limit) const = 0;
    };
    
    
    
} // namespace cyclone

#endif // CONTACTS_H





// SVET CASTIC
#ifndef CYCLONE_PWORLD_H
#define CYCLONE_PWORLD_H

namespace cyclone {
    
    class ParticleWorld
    {
    public:
        typedef std::vector<Particle*> Particles;
        typedef std::vector<ParticleContactGenerator*> ContactGenerators;
        
    protected:

        Particles particles;
        
        bool calculateIterations;
        
        ParticleForceRegistry registry;

        ParticleContactResolver resolver;

        ContactGenerators contactGenerators;

        ParticleContact *contacts;

        unsigned maxContacts;
        
    public:
        
        ParticleWorld(unsigned maxContacts, unsigned iterations=0);

        ~ParticleWorld();
        
        unsigned generateContacts();
        
        void integrate(real duration);
        
        void runPhysics(real duration);
        
        void startFrame();
        
        Particles& getParticles();
        
        ContactGenerators& getContactGenerators();
        
        ParticleForceRegistry& getForceRegistry();
        
        void setMaxContacts(unsigned maxContacts) {
            this->maxContacts = maxContacts;
        }
    };
    
    /**
     * A contact generator that takes an STL vector of particle pointers and
     * collides them against the ground.
     */
    class GroundContacts : public cyclone::ParticleContactGenerator
    {
        cyclone::ParticleWorld::Particles *particles;
        
    public:
        void init(cyclone::ParticleWorld::Particles *particles);
        
        virtual unsigned addContact(cyclone::ParticleContact *contact,
                                    unsigned limit) const;
    };
    
} // namespace cyclone

#endif // PARTICLE_WORLD_H


//################### IMPLEMENTACIA_CORE #######################

using namespace cyclone;

const Vector3 Vector3::GRAVITY = Vector3(0, -9.81, 0);
const Vector3 Vector3::HIGH_GRAVITY = Vector3(0, -19.62, 0);
const Vector3 Vector3::UP = Vector3(0, 1, 0);
const Vector3 Vector3::RIGHT = Vector3(1, 0, 0);
const Vector3 Vector3::OUT_OF_SCREEN = Vector3(0, 0, 1);
const Vector3 Vector3::X = Vector3(1, 0, 0);
const Vector3 Vector3::Y = Vector3(0, 1, 0);
const Vector3 Vector3::Z = Vector3(0, 0, 1);


//################### IMPLEMENTACIA_CORE #######################

//################### IMPLEMENTACIA_PARTICLE #######################

void Particle::integrate(real duration)
{
    // We don't integrate things with zero mass.
    if (inverseMass <= 0.0f) return;
    
    assert(duration > 0.0);
    
    // Update linear position.
    position.addScaledVector(velocity, duration);
    
    // Work out the acceleration from the force
    Vector3 resultingAcc = acceleration;
    resultingAcc.addScaledVector(forceAccum, inverseMass);
    
    // Update linear velocity from the acceleration.
    velocity.addScaledVector(resultingAcc, duration);
    
    // Impose drag.
    velocity *= real_pow(damping, duration);
    
    // Clear the forces.
    clearAccumulator();
}



void Particle::setMass(const real mass)
{
    assert(mass != 0);
    Particle::inverseMass = ((real)1.0)/mass;
}

real Particle::getMass() const
{
    if (inverseMass == 0) {
        return REAL_MAX;
    } else {
        return ((real)1.0)/inverseMass;
    }
}

void Particle::setInverseMass(const real inverseMass)
{
    Particle::inverseMass = inverseMass;
}

real Particle::getInverseMass() const
{
    return inverseMass;
}

bool Particle::hasFiniteMass() const
{
    return inverseMass >= 0.0f;
}

void Particle::setDamping(const real damping)
{
    Particle::damping = damping;
}

real Particle::getDamping() const
{
    return damping;
}

void Particle::setPosition(const Vector3 &position)
{
    Particle::position = position;
}

void Particle::setPosition(const real x, const real y, const real z)
{
    position.x = x;
    position.y = y;
    position.z = z;
}

void Particle::getPosition(Vector3 *position) const
{
    *position = Particle::position;
}

Vector3 Particle::getPosition() const
{
    return position;
}

void Particle::setVelocity(const Vector3 &velocity)
{
    Particle::velocity = velocity;
}

void Particle::setVelocity(const real x, const real y, const real z)
{
    velocity.x = x;
    velocity.y = y;
    velocity.z = z;
}

void Particle::getVelocity(Vector3 *velocity) const
{
    *velocity = Particle::velocity;
}

Vector3 Particle::getVelocity() const
{
    return velocity;
}

void Particle::setAcceleration(const Vector3 &acceleration)
{
    Particle::acceleration = acceleration;
}

void Particle::setAcceleration(const real x, const real y, const real z)
{
    acceleration.x = x;
    acceleration.y = y;
    acceleration.z = z;
}

void Particle::getAcceleration(Vector3 *acceleration) const
{
    *acceleration = Particle::acceleration;
}

Vector3 Particle::getAcceleration() const
{
    return acceleration;
}

void Particle::clearAccumulator()
{
    forceAccum.clear();
}

void Particle::addForce(const Vector3 &force)
{
    forceAccum += force;
}

//################### IMPLEMENTACIA_PARTICLE #######################

//################### IMPLEMENTACIA_PARTICLE_FORCE_GENERATORS #######################

void ParticleForceRegistry::updateForces(real duration)
{
    Registry::iterator i = registrations.begin();
    for (; i != registrations.end(); i++)
    {
        i->fg->updateForce(i->particle, duration);
    }
}

void ParticleForceRegistry::add(Particle* particle, ParticleForceGenerator *fg)
{
    ParticleForceRegistry::ParticleForceRegistration registration;
    registration.particle = particle;
    registration.fg = fg;
    registrations.push_back(registration);
}


ParticleAnchoredSpring::ParticleAnchoredSpring()
{
}

ParticleAnchoredSpring::ParticleAnchoredSpring(Vector3 *anchor,
                                               real sc, real rl)
: anchor(anchor), springConstant(sc), restLength(rl)
{
}

void ParticleAnchoredSpring::init(Vector3 *anchor, real springConstant,
                                  real restLength)
{
    ParticleAnchoredSpring::anchor = anchor;
    ParticleAnchoredSpring::springConstant = springConstant;
    ParticleAnchoredSpring::restLength = restLength;
}



void ParticleAnchoredSpring::updateForce(Particle* particle, real duration)
{
    // Calculate the vector of the spring
    Vector3 force;
    particle->getPosition(&force);
    force -= *anchor;
    
    // Calculate the magnitude of the force
    real magnitude = force.magnitude();
    magnitude = (restLength - magnitude) * springConstant;
    
    // Calculate the final force and apply it
    force.normalise();
    force *= magnitude;
    particle->addForce(force);
}

//################### IMPLEMENTACIA_PARTICLE_FORCE_GENERATORS #######################

//################### IMPLEMENTACIA_PARTICLE_CONTACTS #######################

// Contact implementation

void ParticleContact::resolve(real duration)
{
    resolveVelocity(duration);
    resolveInterpenetration(duration);
}

real ParticleContact::calculateSeparatingVelocity() const
{
    // rychlost, ktorou idu od seba
    Vector3 relativeVelocity = particle[0]->getVelocity();
    if (particle[1]) relativeVelocity -= particle[1]->getVelocity();
    return relativeVelocity * contactNormal;
}

void ParticleContact::resolveVelocity(real duration)
{
    // Find the velocity in the direction of the contact
    real separatingVelocity = calculateSeparatingVelocity();
    
    // Check if it needs to be resolved
    // idu od seba
    if (separatingVelocity > 0)
    {
        // The contact is either separating, or stationary - there's
        // no impulse required.
        return;
    }
    
    // Calculate the new separating velocity
    real newSepVelocity = -separatingVelocity * restitution;
    
    // Check the velocity build-up due to acceleration only
    // rozdiel zrychlenia a nasobene normalou, aby sme vedeli, v ktorom smere to bolo viac
    // trvania - duration
    Vector3 accCausedVelocity = particle[0]->getAcceleration();
    if (particle[1]) accCausedVelocity -= particle[1]->getAcceleration();
    real accCausedSepVelocity = accCausedVelocity * contactNormal * duration;
    
    // If we've got a closing velocity due to acceleration build-up,
    // remove it from the new separating velocity
    if (accCausedSepVelocity < 0)
    {
        newSepVelocity += restitution * accCausedSepVelocity;
        
        // Make sure we haven't removed more than was
        // there to remove.
        if (newSepVelocity < 0) newSepVelocity = 0;
    }
    
    real deltaVelocity = newSepVelocity - separatingVelocity;
    
    // We apply the change in velocity to each object in proportion to
    // their inverse mass (i.e. those with lower inverse mass [higher
    // actual mass] get less change in velocity)..
    real totalInverseMass = particle[0]->getInverseMass();
    if (particle[1]) totalInverseMass += particle[1]->getInverseMass();
    
    // If all particles have infinite mass, then impulses have no effect
    if (totalInverseMass <= 0) return;
    
    // Calculate the impulse to apply
    real impulse = deltaVelocity / totalInverseMass;
    
    // Find the amount of impulse per unit of inverse mass
    Vector3 impulsePerIMass = contactNormal * impulse;
    
    // Apply impulses: they are applied in the direction of the contact,
    // and are proportional to the inverse mass.
    particle[0]->setVelocity(particle[0]->getVelocity() +
                             impulsePerIMass * particle[0]->getInverseMass()
                             );
    if (particle[1])
    {
        // Particle 1 goes in the opposite direction
        particle[1]->setVelocity(particle[1]->getVelocity() +
                                 impulsePerIMass * -particle[1]->getInverseMass()
                                 );
    }
}

void ParticleContact::resolveInterpenetration(real duration)
{
    // If we don't have any penetration, skip this step.
    if (penetration <= 0) return;
    
    // The movement of each object is based on their inverse mass, so
    // total that.
    real totalInverseMass = particle[0]->getInverseMass();
    if (particle[1]) totalInverseMass += particle[1]->getInverseMass();
    
    // If all particles have infinite mass, then we do nothing
    if (totalInverseMass <= 0) return;
    
    // Find the amount of penetration resolution per unit of inverse mass
    Vector3 movePerIMass = contactNormal * (penetration / totalInverseMass);
    
    // Calculate the the movement amounts
    particleMovement[0] = movePerIMass * particle[0]->getInverseMass();
    if (particle[1]) {
        particleMovement[1] = movePerIMass * -particle[1]->getInverseMass();
    } else {
        particleMovement[1].clear();
    }
    
    // Apply the penetration resolution
    particle[0]->setPosition(particle[0]->getPosition() + particleMovement[0]);
    if (particle[1]) {
        particle[1]->setPosition(particle[1]->getPosition() + particleMovement[1]);
    }
}

ParticleContactResolver::ParticleContactResolver(unsigned iterations)
:
iterations(iterations)
{
}

void ParticleContactResolver::setIterations(unsigned iterations)
{
    ParticleContactResolver::iterations = iterations;
}

void ParticleContactResolver::resolveContacts(ParticleContact *contactArray,
                                              unsigned numContacts,
                                              real duration)
{
    unsigned i;
    
    iterationsUsed = 0;
    while(iterationsUsed < iterations)
    {
        // Find the contact with the largest closing velocity;
        real max = REAL_MAX;
        unsigned maxIndex = numContacts;
        for (i = 0; i < numContacts; i++)
        {
            real sepVel = contactArray[i].calculateSeparatingVelocity();
            if (sepVel < max &&
                (sepVel < 0 || contactArray[i].penetration > 0))
            {
                max = sepVel;
                maxIndex = i;
            }
        }
        
        // Do we have anything worth resolving?
        if (maxIndex == numContacts) break;
        
        // Resolve this contact
        contactArray[maxIndex].resolve(duration);
        
        // Update the interpenetrations for all particles
        Vector3 *move = contactArray[maxIndex].particleMovement;
        for (i = 0; i < numContacts; i++)
        {
            if (contactArray[i].particle[0] == contactArray[maxIndex].particle[0])
            {
                contactArray[i].penetration -= move[0] * contactArray[i].contactNormal;
            }
            else if (contactArray[i].particle[0] == contactArray[maxIndex].particle[1])
            {
                contactArray[i].penetration -= move[1] * contactArray[i].contactNormal;
            }
            if (contactArray[i].particle[1])
            {
                if (contactArray[i].particle[1] == contactArray[maxIndex].particle[0])
                {
                    contactArray[i].penetration += move[0] * contactArray[i].contactNormal;
                }
                else if (contactArray[i].particle[1] == contactArray[maxIndex].particle[1])
                {
                    contactArray[i].penetration += move[1] * contactArray[i].contactNormal;
                }
            }
        }
        
        iterationsUsed++;
    }
}

//################### IMPLEMENTACIA_PARTICLE_CONTACTS #######################


//################### IMPLEMENTACIA_PARTICLE_WORLD #######################

ParticleWorld::ParticleWorld(unsigned maxContacts, unsigned iterations)
:
resolver(iterations),
maxContacts(maxContacts)
{
    contacts = new ParticleContact[maxContacts];
    calculateIterations = (iterations == 0);
    
}

ParticleWorld::~ParticleWorld()
{
    delete[] contacts;
}

void ParticleWorld::startFrame()
{
    for (Particles::iterator p = particles.begin();
         p != particles.end();
         p++)
    {
        // Remove all forces from the accumulator
        (*p)->clearAccumulator();
    }
}

unsigned ParticleWorld::generateContacts()
{
    unsigned limit = maxContacts;
    ParticleContact *nextContact = contacts;
    
    for (ContactGenerators::iterator g = contactGenerators.begin();
         g != contactGenerators.end();
         g++)
    {
        unsigned used =(*g)->addContact(nextContact, limit);
        limit -= used;
        nextContact += used;
        
        // We've run out of contacts to fill. This means we're missing
        // contacts.
        if (limit <= 0) break;
    }
    
    // Return the number of contacts used.
    return maxContacts - limit;
}

void ParticleWorld::integrate(real duration)
{
    for (Particles::iterator p = particles.begin();
         p != particles.end();
         p++)
    {
        // Remove all forces from the accumulator
        (*p)->integrate(duration);
    }
}

void ParticleWorld::runPhysics(real duration)
{
    // First apply the force generators
    registry.updateForces(duration);
    
    // Then integrate the objects
    integrate(duration);
    
    // Generate contacts
    unsigned usedContacts = generateContacts();
    
    // And process them
    if (usedContacts)
    {
        if (calculateIterations) resolver.setIterations(usedContacts * 2);
        resolver.resolveContacts(contacts, usedContacts, duration);
    }
}


ParticleWorld::Particles& ParticleWorld::getParticles()
{
    return particles;
}

ParticleWorld::ContactGenerators& ParticleWorld::getContactGenerators()
{
    return contactGenerators;
}

ParticleForceRegistry& ParticleWorld::getForceRegistry()
{
    return registry;
}

void GroundContacts::init(cyclone::ParticleWorld::Particles *particles)
{
    GroundContacts::particles = particles;
}

unsigned GroundContacts::addContact(cyclone::ParticleContact *contact,
                                    unsigned limit) const
{
    unsigned count = 0;
    for (cyclone::ParticleWorld::Particles::iterator p = particles->begin();
         p != particles->end();
         p++)
    {
        cyclone::real y = (*p)->getPosition().y;
        if (y < 0.01f+(*p)->getSize())
        {
            contact->contactNormal = cyclone::Vector3::UP;
            contact->particle[0] = *p;
            contact->particle[1] = NULL;
            contact->penetration = -y;
            contact->restitution = 0.4f;
            contact++;
            count++;
        }
        
        if (count >= limit) return count;
    }
    return count;
}
//################### IMPLEMENTACIA_PARTICLE_WORLD #######################


// APPLICATION

class Application
{
protected:

public:

    int PLAYGROUND_HEIGHT;

    int PLAYGROUND_WIDTH;
    
    virtual const char* getTitle();
    
    virtual void initGraphics();

    virtual void setView();

    virtual void deinit();

    virtual void display();

    virtual void update();

    virtual void key(int key);
    
    virtual void keyUP(int key);
    
    virtual void resize(int width, int height);
    
};

class MassAggregateApplication : public Application
{
protected:
    cyclone::ParticleWorld world;
    cyclone::Particle *particleArray;
    cyclone::GroundContacts groundContactGenerator;
    
public:
    MassAggregateApplication(unsigned int particleCount);
    virtual ~MassAggregateApplication();
    
    virtual void update();
    
    virtual void initGraphics();
    
    virtual void display();
};
// APPLICATION_HEADER_END


// TIMING
#ifndef CYCLONE_DEMO_TIMING_H
#define CYCLONE_DEMO_TIMING_H

struct TimingData
{
    unsigned frameNumber;
    
    unsigned lastFrameTimestamp;

    unsigned lastFrameDuration;

    unsigned long lastFrameClockstamp;
    
    unsigned long lastFrameClockTicks;
    
    bool isPaused;

    double averageFrameDuration;
    
    float fps;
    
    static TimingData& get();

    static void update();

    static void init();
    
    static void deinit();
    
    static unsigned getTime();
    
    static unsigned long getClock();
    
    
private:
    // These are private to stop instances being created: use get().
    TimingData() {}
    TimingData(const TimingData &) {}
    TimingData& operator=(const TimingData &);
};


#endif // TIMING_H

// ################ IMPLEMENTACIA_TIMING######################
// Hold internal timing data for the performance counter.
static bool qpcFlag;

#if (__APPLE__ || __unix)
#define TIMING_UNIX	1

#include <stdlib.h>
#include <sys/time.h>

// assume unix based OS
typedef unsigned long long	LONGLONG;
#else
#define TIMING_WINDOWS	1
// assume windows

// Import the high performance timer (c. 4ms).
#include <windows.h>
#include <mmsystem.h>

static double qpcFrequency;
#endif



// Internal time and clock access functions
unsigned systemTime()
{
#if TIMING_UNIX
    struct timeval tv;
    gettimeofday(&tv, 0);
    
    return (unsigned)tv.tv_sec * 1000 + tv.tv_usec/1000;
    
#else
    if(qpcFlag)
    {
        static LONGLONG qpcMillisPerTick;
        QueryPerformanceCounter((LARGE_INTEGER*)&qpcMillisPerTick);
        return (unsigned)(qpcMillisPerTick * qpcFrequency);
    }
    else
    {
        return unsigned(timeGetTime());
    }
#endif
    
}

unsigned TimingData::getTime()
{
    return systemTime();
}

#if TIMING_WINDOWS
unsigned long systemClock()
{
    __asm {
        rdtsc;
    }
}
#endif

unsigned long TimingData::getClock()
{
    
#if TIMING_UNIX
    struct timeval tv;
    gettimeofday(&tv, 0);
    
    return tv.tv_sec * 1000 + tv.tv_usec/1000;
#else
    return systemClock();
#endif
}

// Sets up the timing system and registers the performance timer.
void initTime()
{
#if TIMING_UNIX
    qpcFlag = false;
#else
    LONGLONG time;
    
    qpcFlag = (QueryPerformanceFrequency((LARGE_INTEGER*)&time) > 0);
    
    // Check if we have access to the performance counter at this
    // resolution.
    if (qpcFlag) qpcFrequency = 1000.0 / time;
#endif
}


// Holds the global frame time that is passed around
static TimingData *timingData = NULL;

// Retrieves the global frame info instance
TimingData& TimingData::get()
{
    return (TimingData&)*timingData;
}

// Updates the global frame information. Should be called once per frame.
void TimingData::update()
{
    if (!timingData) return;
    
    // Advance the frame number.
    if (!timingData->isPaused)
    {
        timingData->frameNumber++;
    }
    
    // Update the timing information.
    unsigned thisTime = systemTime();
    timingData->lastFrameDuration = thisTime -
    timingData->lastFrameTimestamp;
    timingData->lastFrameTimestamp = thisTime;
    
    // Update the tick information.
    unsigned long thisClock = getClock();
    timingData->lastFrameClockTicks =
    thisClock - timingData->lastFrameClockstamp;
    timingData->lastFrameClockstamp = thisClock;
    
    // Update the RWA frame rate if we are able to.
    if (timingData->frameNumber > 1) {
        if (timingData->averageFrameDuration <= 0)
        {
            timingData->averageFrameDuration =
            (double)timingData->lastFrameDuration;
        }
        else
        {
            // RWA over 100 frames.
            timingData->averageFrameDuration *= 0.99;
            timingData->averageFrameDuration +=
            0.01 * (double)timingData->lastFrameDuration;
            
            // Invert to get FPS
            timingData->fps =
            (float)(1000.0/timingData->averageFrameDuration);
        }
    }
}

void TimingData::init()
{
    // Set up the timing system.
    initTime();
    
    // Create the frame info object
    if (!timingData) timingData = new TimingData();
    
    // Set up the frame info structure.
    timingData->frameNumber = 0;
    
    timingData->lastFrameTimestamp = systemTime();
    timingData->lastFrameDuration = 0;
    
    timingData->lastFrameClockstamp = getClock();
    timingData->lastFrameClockTicks = 0;
    
    timingData->isPaused = false;
    
    timingData->averageFrameDuration = 0;
    timingData->fps = 0;
}

void TimingData::deinit()
{
    delete timingData;
    timingData = NULL;
}
// ################ IMPLEMENTACIA_TIMING######################

// ################ IMPLEMENTACIA_APPLICATION######################

void Application::initGraphics()
{
    glClearColor(0.9f, 0.95f, 1.0f, 1.0f);
    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_SMOOTH);
    
    setView();
}

void Application::setView()
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60.0, (double)PLAYGROUND_WIDTH/(double)PLAYGROUND_HEIGHT, 1.0, 500.0);
    glMatrixMode(GL_MODELVIEW);
}

void Application::display()
{
    glClear(GL_COLOR_BUFFER_BIT);
    
    glBegin(GL_LINES);
    glVertex2i(1, 1);
    glVertex2i(639, 319);
    glEnd();
}

const char* Application::getTitle()
{
    return "Defaultne meno";
}

void Application::deinit()
{
}

void Application::update()
{
    glutPostRedisplay();
}

void Application::key(int key)
{
}

void Application::keyUP(int key) {
    
}


void Application::resize(int width, int height)
{
    // Avoid the divide by zero.
    if (height <= 0) height = 1;
    
    // Set the internal variables and update the view
    Application::PLAYGROUND_WIDTH = width;
    Application::PLAYGROUND_HEIGHT = height;
    glViewport(0, 0, width, height);
    setView();
}


MassAggregateApplication::MassAggregateApplication(unsigned int particleCount)
:
world(particleCount*10)
{
    particleArray = new cyclone::Particle[particleCount];
    for (unsigned i = 0; i < particleCount; i++)
    {
        world.getParticles().push_back(particleArray + i);
    }
    
    groundContactGenerator.init(&world.getParticles());
    world.getContactGenerators().push_back(&groundContactGenerator);
}

MassAggregateApplication::~MassAggregateApplication()
{
    delete[] particleArray;
}

void MassAggregateApplication::initGraphics()
{
    // Call the superclass
    Application::initGraphics();
}

void MassAggregateApplication::display()
{

}

void MassAggregateApplication::update()
{
    // Clear accumulators
    world.startFrame();
    
    // Find the duration of the last frame in seconds
    float duration = (float)TimingData::get().lastFrameDuration * 0.001f;
    if (duration <= 0.0f) return;
    
    // Run the simulation
    world.runPhysics(duration);
    
    Application::update();
}

// ################ IMPLEMENTACIA_APPLICATION######################


//################## MAIN ##################################################

extern Application* getApplication();

// Store the global application object.
Application* app;

void createWindow(const char* title)
{
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(1280,640);
    glutInitWindowPosition(0,0);
    glutCreateWindow(title);
}

void update()
{
    // Update the timing.
    TimingData::get().update();
    
    // Delegate to the application.
    app->update();
}

void display()
{
    app->display();
    
    // Update the displayed content.
    glFlush();
    glutSwapBuffers();
}

void reshape(int width, int height)
{
    app->resize(width, height);
}

void keyboard(int key, int x, int y)
{
    // Note we omit passing on the x and y: they are rarely needed.
    app->key(key);
}

void keyboardUP(int key, int x, int y) {
    app->keyUP(key);
}

int main(int argc, char** argv)
{
    // Set up GLUT and the timers
    glutInit(&argc, argv);
    TimingData::init();
    
    // Create the application and its window
    app = getApplication();
    createWindow(app->getTitle());
    
    // Set up the appropriate handler functions
    glutReshapeFunc(reshape);
    glutSpecialFunc(keyboard);
    glutSpecialUpFunc(keyboardUP);
    glutDisplayFunc(display);
    glutIdleFunc(update);
    

    app->initGraphics();
    glutMainLoop();
    
    app->deinit();
    delete app;
    TimingData::deinit();
}
//################## MAIN ####################


//##################### SIMULATION CODE ############################

//## GLOBALS ##
#define GAME_CHAR_INDEX 0
#define GAMECHAR_SIZE 0.2f
#define CAMERA_INDEX 1
#define SPACE_KEY 32
#define DEFAULT_PARTICLES 2
#define PLAYGROUND_HEIGHT 64
#define PLAYGROUND_WIDTH 64
#define BASE_MASS 1
#define BASE_DAMPING 0.75f

//## COLORS ###
const Vector3 white = Vector3(1,1,1);
const Vector3 black = Vector3(0,0,0);
const Vector3 green = Vector3(0,1,0);
const Vector3 red = Vector3(1,0,0);
const GLfloat shadow_gray[] = {0.0, 0.0, 0.0, 0.8};
const float DEG2RAD = 3.14159/180;
bool playingGame = false;
#define WHITE (GLfloat *)&white
#define BLACK (GLfloat *)&black
#define RED (GLfloat *)&red
#define GREEN (GLfloat *)&green
#define SHADOW (GLfloat *)&shadow_gray


std::string convertInt(int number) {
    std::stringstream ss;//create a stringstream
    ss << number;//add number to the stream
    return ss.str();//return a string with the contents of the stream
}

struct GameObject{
    Vector3 pos;
    real size;
};

class Draw {
public:
    static void square(GLfloat x, GLfloat y, GLfloat z, GLfloat size = 1) {
        glBegin(GL_QUADS);
        glNormal3f(0.0f,-1.0f,0.0f);
        glVertex3f(x, y, z);
        glNormal3f(0.0f,-1.0f,0.0f);
        glVertex3f(x+size, y, z);
        glNormal3f(0.0f,-1.0f,0.0f);
        glVertex3f(x+size, y, z+size);
        glNormal3f(0.0f,-1.0f,0.0f);
        glVertex3f(x, y, z+size);
        glEnd();
    }
    
    static void circle(float radius)
    {
        glBegin(GL_POLYGON);
        
        for (int i=0; i < 360; i++)
        {
            float degInRad = i*DEG2RAD;
            glVertex3f(cos(degInRad)*radius,0.0f,sin(degInRad)*radius);
        }
        
        glEnd();
    }
    static void renderBitmapString(float x, float y, void *font, std::string s) {
        
        //TEXT
        glMatrixMode( GL_PROJECTION ) ;
        glPushMatrix() ;
        glLoadIdentity();
        glMatrixMode( GL_MODELVIEW ) ;
        glPushMatrix() ;
        glLoadIdentity() ;
        glDisable(GL_LIGHTING);
        
        glRasterPos2f( x,y ) ;
        
        for (int i=0; i < s.size(); i++)
        {
            glutBitmapCharacter(font, s[i]);
        }
        
        glEnable(GL_LIGHTING);
        glMatrixMode( GL_PROJECTION ) ;
        glPopMatrix() ;
        glMatrixMode( GL_MODELVIEW ) ;
        glPopMatrix() ;
    }
};


class WallContacts : public ParticleContactGenerator {
// nastav kontakty so stenami ihriska
protected:
    ParticleWorld::Particles *particles;
public:
    void init(ParticleWorld::Particles *particles) {
        this->particles = particles;
    }
    
    unsigned addContact(ParticleContact *contact, unsigned limit) const {
        
        
        unsigned count = 0;
        for (cyclone::ParticleWorld::Particles::iterator p = particles->begin();
             p != particles->end();
             p++)
        {
            real x = (*p)->getPosition().x;
            real z = (*p)->getPosition().z;
            if (x < -PLAYGROUND_WIDTH/2 || x > PLAYGROUND_WIDTH/2) {
                contact->contactNormal = (x < 0) ? Vector3::X : Vector3::X*(-1);
                contact->particle[0] = *p;
                contact->particle[1] = NULL;
                contact->penetration = (x < 0) ? x + PLAYGROUND_WIDTH/2 : x - PLAYGROUND_WIDTH/2;
                contact->restitution = 0.8f;
                contact++;
                count++;
            } else if (z < -PLAYGROUND_HEIGHT/2 || z > PLAYGROUND_HEIGHT/2) {
                contact->contactNormal = (z < 0) ? Vector3::Z : Vector3::Z*(-1);
                contact->particle[0] = *p;
                contact->particle[1] = NULL;
                contact->penetration = (z < 0) ? z + PLAYGROUND_HEIGHT/2 : z - PLAYGROUND_HEIGHT/2;
                contact->restitution = 0.8f;
                contact++;
                count++;
            }
            if (count >= limit) return count;
        }
        return count;
    }
    
};

class EnemyContacts : public ParticleContactGenerator {
protected:
    ParticleWorld::Particles *particles;
public:
    
    void init(ParticleWorld::Particles *particles) {
        this->particles = particles;
    }
    
    unsigned addContact(ParticleContact *contact, unsigned limit) const {
        unsigned count = 0;
        cyclone::ParticleWorld::Particles::iterator p1 = particles->begin();
        cyclone::ParticleWorld::Particles::iterator p2 = particles->begin();
        for (; p1 != particles->end(); p1++) {
            // dumb check ci je to gamechar particle - lebo tam chceme hned koniec hry
            if ((*p1)->getSize() == GAMECHAR_SIZE && playingGame) {
                continue;
            }
            real p1Size = (*p1)->getSize();
            for (; p2 != particles->end(); p2++) {
                // dumb check ci je to gamechar particle - lebo tam chceme hned koniec hry
                if (((*p1) == (*p2) || (*p2)->getSize() == GAMECHAR_SIZE) && playingGame) {
                    continue;
                }
                
                real distBetween = ((*p1)->getPosition() - (*p2)->getPosition()).magnitude();
                if (distBetween < p1Size + (*p2)->getSize()) {
                    Vector3 contNormal = (*p1)->getPosition() - (*p2)->getPosition();
                    contNormal.normalise();
                    
                    contact->contactNormal = contNormal;
                    contact->particle[0] = *p1;
                    contact->particle[1] = *p2;
                    contact->penetration = (p1Size + (*p2)->getSize()) - distBetween;
                    contact->restitution = 0.8f;
                    contact++;
                    count++;
                }
            }
            p2 = particles->begin();
            
            
            if (count >= limit) return count;
        }
        
        return count;
    }
    
};

class CameraSpring : public ParticleAnchoredSpring {
// pruzina, ktora sposobuje pohyb kamery
public:
    CameraSpring(Vector3 *anchor, real springConst, real restLength):ParticleAnchoredSpring(anchor, springConst, restLength){}
    void changeAnchorPos(Vector3 newAnchor) {
        *anchor = newAnchor;
    }
    void updateForce(Particle* particle, real duration)
    {
        Vector3 force;
        particle->getPosition(&force);
        force -= *anchor;
        
        real magnitude = force.magnitude();
        if (magnitude < restLength) return;
        
        magnitude = magnitude - restLength;
        magnitude *= springConstant;
        
        force.normalise();
        force *= -magnitude;
        
        // ak by mala zacat kamera klesat - uplatnime silu v smere y
        if (anchor->y + 2.0f > particle->getPosition().y) {
            force.y = 4.0f;
        }
        particle->addForce(force);
    }
};
class Camera {
public:
    Particle *particleRef;
    CameraSpring *cameraLink;
public:
    Camera(Particle *p): particleRef(p) {
        particleRef->setAcceleration(0,0,0);
        particleRef->setDamping(0.1f);
        particleRef->setMass(1);
        particleRef->setPosition(5,1.5,10.0);
    }
    void respawn() {
        particleRef->setPosition(30,10.5,10.0);
    }
};
class GameCharacter {
/**
 * Nas hlavny herny charakter
 * spravanie : particle
 * damping : zmenseny an 0.75f - rozumna hodnota, kedy sa castica nechova ako vo vode
 * + budi dojem spravneho trenia
 **/
public:
    Particle *particleRef;
    Vector3 *charPos;
    const real charForward = 8;
    const real charBackward = 8;
    const real charSteering = 12;
    real size = GAMECHAR_SIZE;
    int lives = 3;
    bool isInvincible = false;
    unsigned hitTime = -1;
    bool isHit = false;
protected:

    Vector3 movementAccel = Vector3::GRAVITY;
    // constants
    const real jumpForce = 6;
public:
    GameCharacter(Particle *particle): particleRef(particle) {
        particleRef->setMass(BASE_MASS);
        particleRef->setAcceleration(Vector3::GRAVITY);
        particleRef->setDamping(BASE_DAMPING);
        particleRef->setPosition(0,2,5);
        particleRef->setSize(size);
        charPos = new Vector3(particleRef->getPosition());
    }
    
    /**
     * vykreslovanie charakteru
     **/
    void draw() {
        Vector3 pos = getPosition();
        glPushMatrix();
            glColor3fv(GREEN);
            if (isHit) {
                if (((TimingData::getTime() - hitTime)*0.001f) < 2.00f) {
                    if(fmod(((TimingData::getTime() - hitTime)*0.001f), 0.2) < 0.1) {
                        glColor3f(1.0, 0.5, 1.0);
                    }
                } else {
                    this->hitTime = -1.0f;
                    this->isInvincible = false;
                    this->isHit = false;
                }
            }
        
            glTranslatef(pos.x, pos.y, pos.z);
            glutSolidSphere(size,20,20);
            glTranslatef(0.0f, -pos.y+0.01f, 0.0);
            glColor4fv(SHADOW);
            Draw::circle(size);
        glPopMatrix();
        //shadow

        
    }
    
    
    
    
    /**
     * Pohyb postavy
     * zavisi od kamery - ukazalo sa, ze pre pohyb lopticky, je toto
     * zatial najlepsi sposob pohybu - kedze je kamera ukotvena ako pruzina
     * pohyb dopredu je vzdy v smere kamera lopticka
     **/
    void move(map<int, real> keysPressed, Camera *camera) {
        Vector3 vel = particleRef->getVelocity();
        vel.normalise();
        Vector3 force = Vector3::GRAVITY;
        map<int,real>::iterator j = keysPressed.begin();
        for (; j != keysPressed.end(); j++) {
            switch (j->first) {
                case 'w': {
                    if (j->second != 0.0f) {
                        Vector3 cameraCharNormal = camera->particleRef->getPosition() - particleRef->getPosition();
                        cameraCharNormal.normalise();
                        cameraCharNormal.y =0;
                        cameraCharNormal.invert();
                        Vector3 forwardVel = particleRef->getVelocity();
                        forwardVel.addScaledVector(cameraCharNormal,0.15);
                        particleRef->setVelocity(forwardVel);
                    }
                }
                    break;
                case 's': {
                    if (j->second != 0.0f) {
                        Vector3 cameraCharNormal = camera->particleRef->getPosition() - particleRef->getPosition();
                        cameraCharNormal.normalise();
                        cameraCharNormal.y =0;
                        Vector3 backwardVel = particleRef->getVelocity();
                        backwardVel.addScaledVector(cameraCharNormal,0.2);
                        particleRef->setVelocity(backwardVel);
                    }
                }
                    break;
                case 'a': {
                    if (j->second != 0.0f) {
                        Vector3 cameraCharNormal = camera->particleRef->getPosition() - particleRef->getPosition();
                        cameraCharNormal.normalise();
                        Vector3 steerLeft = cameraCharNormal.vectorProduct(Vector3::UP);
                        real newx,newz;
                        newx = steerLeft.x * cos(-M_PI/4) - steerLeft.z * sin(-M_PI/4);
                        newz = steerLeft.x * sin(-M_PI/4) + steerLeft.z * cos(-M_PI/4);
                        steerLeft.x = newx;
                        steerLeft.z = newz;
                        Vector3 steerVel = particleRef->getVelocity();
                        steerVel.addScaledVector(steerLeft,0.3);
                        particleRef->setVelocity(steerVel);
                    }
                }
                    break;
                case 'd':
                    if (j->second != 0.0f) {
                        Vector3 cameraCharNormal = camera->particleRef->getPosition() - particleRef->getPosition();
                        cameraCharNormal.normalise();
                        Vector3 steerRight = Vector3::UP.vectorProduct(cameraCharNormal);
                        real newx,newz;
                        newx = steerRight.x * cos(M_PI/4) - steerRight.z * sin(M_PI/4);
                        newz = steerRight.x * sin(M_PI/4) + steerRight.z * cos(M_PI/4);
                        steerRight.x = newx;
                        steerRight.z = newz;
                        Vector3 steerVel = particleRef->getVelocity();
                        steerVel.addScaledVector(steerRight,0.3);
                        particleRef->setVelocity(steerVel);
                    }
                    break;
            }
        }
        particleRef->setAcceleration(force);
    }

    void jump() {
        if (particleRef->getPosition().y < particleRef->getSize() + 0.2f) {
            Vector3 newVel = particleRef->getVelocity() + Vector3::UP*jumpForce;
            particleRef->setVelocity(newVel);
        }
    }
    
    void update() {
        charPos = new Vector3(particleRef->getPosition());
    }
    
    void respawn() {
        particleRef->setVelocity(0, 0, 0);
        particleRef->setPosition(0,10,5);
        particleRef->setAcceleration(Vector3::GRAVITY);
        lives = 3;
        hitTime = -1.0f;
        isInvincible = false;
        isHit = false;
    }

    Vector3 getPosition() {
        return particleRef->getPosition();
    }
    real getSize() {
        return particleRef->getSize();
    }
    
    void wasHit(unsigned hitTime) {
        this->hitTime = hitTime;
        this->isInvincible = true;
        this->isHit = true;
    }

};


class GamePoint{
protected:
    GLfloat x;
    GLfloat y;
    GLfloat z;
    int size;
    bool acquired;
    bool visible;
    int score;
    unsigned timeAcquired;
public:
    GamePoint() {
        x = rand() % ((PLAYGROUND_WIDTH/2)-1) * (rand()%2 ? -1 : 1);
        z = rand() % ((PLAYGROUND_HEIGHT/2)-1) * (rand()%2 ? -1 : 1);
        y = 0.001f;
        size = 2;
        acquired = false;
        visible = true;
        score = 1;
    };
    void drawPoint() {
        if (visible){
            if ((acquired && ((TimingData::getTime() - timeAcquired)*0.001f) > 2.00f)) {
                visible = false;
                return;
            }
            glColor4f(1.0, 0.0, 0.0, 0.8);
            
            if (acquired) {

                if (fmod(((TimingData::getTime() - timeAcquired)*0.001f), 0.6) < 0.3) {
                    glColor4f(0.0, 1.0, 0.0, 0.8);
                }
                
            }
            glPushMatrix();
                Draw::square(x, y, z, size);
            glPopMatrix();
        }
    };
    
    int update(Vector3 charPos) {
        if (!visible) {
            refresh();
        }
        if (acquired) {
            return 0;
        }
        if (charPos.x > x && charPos.x < x+size && charPos.z > z && charPos.z < z+size) {
            acquired = true;
            timeAcquired = TimingData::getTime();
            return score;
        }
        return 0;
    }
    
    bool isAcquired() const{
        return acquired;
    }
    
    void refresh() {
        x = rand() % ((PLAYGROUND_WIDTH/2)-1) * (rand()%2 ? -1 : 1);
        z = rand() % ((PLAYGROUND_HEIGHT/2)-1) * (rand()%2 ? -1 : 1);
        acquired = false;
        visible = true;
    }
};

class Enemy {
protected:
    real x;
    real y;
    real z;
    GLfloat size = 0.1f;
    Vector3 color;
    bool hitPlayer = false;
public:
    Particle *particleRef;
public:
    Enemy(real x, real y, real z, Particle *p): x(x),y(y),z(z),particleRef(p){
        p->setPosition(Vector3(x,y,z));
        p->setAcceleration(Vector3::GRAVITY);
        p->setMass(BASE_MASS);
        p->setDamping(0.995);
        p->setSize(size);
        p->clearAccumulator();
    }
    void draw() {
        glPushMatrix();
        Vector3 pos = getPosition();
        glColor3fv((GLfloat *)&color);
        glTranslatef(pos.x, pos.y, pos.z);
        glutSolidSphere(size, 10, 20);
        glTranslatef(0.0f, -pos.y+0.01f, 0.0);
        glColor4fv(SHADOW);
        Draw::circle(size);
        glPopMatrix();
    }
    virtual void move(GameObject charPos) = 0;
    
    virtual void update(GameObject charPos) {
        Enemy::hitPlayer = false;
        checkHitPlayer(charPos);
        move(charPos);
    }
    Vector3 getPosition() {
        return particleRef->getPosition();
    }
    bool didHitPlayer() {
        return hitPlayer;
    }
    virtual void checkHitPlayer(GameObject charInfo) {
        Vector3 enemyCharVec = getPosition() - charInfo.pos;
        if (enemyCharVec.magnitude() <= charInfo.size + size) {
            Enemy::hitPlayer = true;
        }
    }
};

/**
 * Nepriatel pohybujuci sa po jednej osi
 **/
class EnemyLVL1:public Enemy {
protected:
    Vector3 startPos;
    
    const GLfloat size = 0.5;
    const real allowedDistance = 5.0f;
    const real speed = 10.0f;
    Vector3 color = Vector3(255.0f/255.0f,145.0f/255.0f,171.0f/255.0f);
public:
    EnemyLVL1(real x, real y, real z,Particle *p):Enemy(x,y,z,p){
        startPos = Vector3(x,y,z);
        p->setVelocity(rand()%(int)speed,0.0f,0);
        Enemy::size = size;
        Enemy::color = color;
        p->setSize(size);
    }
    void move(GameObject charPos) {
        Vector3 vel = particleRef->getVelocity();
        vel.normalise();
        vel.addScaledVector(vel, speed);
        particleRef->setVelocity(vel);
    }
    
};

/**
 * Nepriatel pohybujuci sa po osi X a Z zaroven
 **/
class EnemyLVL2:public Enemy {
protected:
    Vector3 startPos;
    
    const GLfloat size = 0.5;
    const real allowedDistance = 5.0f;
    const real speed = 10.0f;
    Vector3 color = Vector3(1,59.0f/255.0f,104.0f/255.0f);
public:
    EnemyLVL2(real x, real y, real z,Particle *p):Enemy(x,y,z,p){
        startPos = Vector3(x,y,z);
        p->setVelocity(rand()%(int)speed,0.0f,rand()%(int)speed);
        Enemy::size = size;
        Enemy::color = color;
        p->setSize(size);
    }
    void move(GameObject charPos) {
        Vector3 vel = particleRef->getVelocity();
        vel.normalise();
        vel.addScaledVector(vel, speed);
        particleRef->setVelocity(vel);
    }
};

class EnemyLVL3:public Enemy {
protected:
    Vector3 startPos;
    
    const GLfloat size = 0.5;
    real changeDirectionDist;
    const real speed = 10.0f;
    unsigned timeAppeared;
    Vector3 color = Vector3(1,0,0);
public:
    EnemyLVL3(real x, real y, real z,Particle *p):Enemy(x,y,z,p){
        startPos = Vector3(x,y,z);
        p->setVelocity(rand()%(int)speed,0.0f,rand()%(int)speed);
        Enemy::size = size;
        Enemy::color = color;
        p->setSize(size);
        
        changeDirectionDist = 10.0f;
    }
    void move(GameObject charPos) {
        Vector3 dir = startPos - getPosition();
        real distance = dir.magnitude();
        if (distance > changeDirectionDist) {
            Vector3 newDir = Vector3((rand() % 5 + 5) * (rand() % 2 ? -1 : 1),0,(rand() % 5 + 5) * (rand() % 2 ? -1 : 1));
            particleRef->setVelocity(newDir);
            startPos = getPosition();
        }
        Vector3 vel = particleRef->getVelocity();
        vel.normalise();
        vel.addScaledVector(vel, speed);
        particleRef->setVelocity(vel);
    }
    
};

class EnemyLVL4:public Enemy {
protected:
    Vector3 startPos;
    
    const GLfloat size = 0.5;
    const real speed = 5.0f;
    unsigned timeAppeared;
    Vector3 color = Vector3(87.0f/255.0f,0,20.0f/255.0f);
public:
    EnemyLVL4(real x, real y, real z,Particle *p):Enemy(x,y,z,p){
        startPos = Vector3(x,y,z);
        Enemy::size = size;
        Enemy::color = color;
        p->setSize(size);
        p->setVelocity(0.0f,0.0f,0.0f);
    }
    void move(GameObject charPos) {
        Vector3 dirToChar = charPos.pos - getPosition();
        dirToChar.normalise();
        dirToChar.y = 0;
        Vector3 vel = particleRef->getVelocity();
        vel.addScaledVector(dirToChar, 0.05f);
        particleRef->setVelocity(vel);
    }
    
};

class Playground {
public:
    unsigned gameDuration;
    int score;
    bool endGame = false;
protected:
    GameCharacter *gamechar;
    Camera *camera;
    ParticleWorld *world;
    vector<GamePoint*> points;
    vector<Enemy*> enemies;
    int width = PLAYGROUND_WIDTH;
    int height = PLAYGROUND_HEIGHT;
    int gamePoints = 5;
    int levelChangeFactor = 5;
    bool inMenu = true;
protected:
    void drawFloor() {
        int i,j;
        for (i = 0; i < height; i++) {
            for (j = 0; j < width; j++) {
                if ((j%2)^(i%2)) {
                    glColor3fv(WHITE);
                } else {
                    glColor3fv(BLACK);
                }
                Draw::square(-(width/2)+j, 0, -(height/2)+i);
            }
        }
    }
    
    void drawPoints() {
        vector<GamePoint*>::iterator it = points.begin();
        for (; it != points.end(); it++) {
            (*it)->drawPoint();
        }
    }
    
    void drawEnemies() {
        vector<Enemy*>::iterator en = enemies.begin();
        for (; en != enemies.end(); en++) {
            (*en)->draw();
        }
    }
    
    void updatePoints() {
        vector<GamePoint*>::iterator it = points.begin();
        int addPoint = 0;
        for (; it != points.end(); it++) {
            if((addPoint = (*it)->update(gamechar->getPosition()))){
                score += addPoint;
                // pred kazdou zmenou levelu vymaze minuly level superov
                if (score % levelChangeFactor == 0) {
                    removeEnemies();
                }
                return;
            }
        }
    }
    void updateEnemies() {
        vector<Enemy *>::iterator en = enemies.begin();
        for (; en != enemies.end(); en++) {
            GameObject gamecharInfo;
            gamecharInfo.pos = gamechar->getPosition();
            gamecharInfo.size = gamechar->getSize();
            (*en)->update(gamecharInfo);
            if ((*en)->didHitPlayer() && !gamechar->isInvincible && playingGame) {
                gamechar->lives--;
                gamechar->wasHit(TimingData::getTime());
            }
        }
        
        
    }
    
private:
    
    /**
     * Vrati nahodnu poziciu na hracej ploche
     **/
    Vector3 getRandomPlaygroundPos() {
        return Vector3(rand() % PLAYGROUND_WIDTH/2 * (rand()%2 ? -1 : 1),0,rand() % PLAYGROUND_HEIGHT/2 * (rand()%2 ? -1 : 1));
    }

    Vector3 getSafePosition() {
        Vector3 enemyPos = getRandomPlaygroundPos();
        real enemyCharDist = (gamechar->getPosition() - getRandomPlaygroundPos()).magnitude();
        while (enemyCharDist < 5.0f) {
            enemyPos = getRandomPlaygroundPos();
            enemyCharDist = (gamechar->getPosition() - getRandomPlaygroundPos()).magnitude();
        }
        
        return enemyPos;
    }
    
    /**
     * Vygeneruje nepriatelov LEVELU 1 - EnemyLVL1
     **/
    

    
    void spawnLVL1Enemies(int quantity) {

        if (enemies.size() == 0) {
            for (int i = 0; i < quantity; i++) {
                world->getParticles().push_back(new Particle());
                Vector3 enemyPos = getSafePosition();
                enemies.push_back(new EnemyLVL1(enemyPos.x,3,enemyPos.z,world->getParticles().back()));
            }
        }
        
    }
    
    void spawnLVL2Enemies(int quantity) {
        if (enemies.size() == 0) {
            for (int i = 0; i < quantity; i++) {
                world->getParticles().push_back(new Particle());
                Vector3 enemyPos = getSafePosition();
                enemies.push_back(new EnemyLVL2(enemyPos.x,3,enemyPos.z,world->getParticles().back()));
            }
        }
    }
    
    void spawnLVL3Enemies(int quantity) {
        if (enemies.size() == 0) {
            for (int i = 0; i < quantity; i++) {
                world->getParticles().push_back(new Particle());
                Vector3 enemyPos = getSafePosition();
                enemies.push_back(new EnemyLVL3(enemyPos.x,3,enemyPos.z,world->getParticles().back()));
            }
        }
    }
    
    void spawnLVL4Enemies(int quantity) {
        if (enemies.size() == 0) {
            for (int i = 0; i < quantity; i++) {
                world->getParticles().push_back(new Particle());
                Vector3 enemyPos = getSafePosition();
                enemies.push_back(new EnemyLVL4(enemyPos.x,3,enemyPos.z,world->getParticles().back()));
            }
            
        }
    }
    void spawnLVL5Enemies(int quantity) {
        if (enemies.size() == 0) {
            for (int i = 0; i < quantity; i++) {
                world->getParticles().push_back(new Particle());
                Vector3 enemyPos = getSafePosition();
                switch (i % 2) {
                    case 0:
                        enemies.push_back(new EnemyLVL1(enemyPos.x,3,enemyPos.z,world->getParticles().back()));
                        break;
                    case 1:
                        enemies.push_back(new EnemyLVL2(enemyPos.x,3,enemyPos.z,world->getParticles().back()));
                        break;
                }

            }
        }
    }
    void spawnLVL6Enemies(int quantity) {
        if (enemies.size() == 0) {
            for (int i = 0; i < quantity; i++) {
                world->getParticles().push_back(new Particle());
                Vector3 enemyPos = getSafePosition();
                switch (i % 2) {
                    case 0:
                        enemies.push_back(new EnemyLVL1(enemyPos.x,3,enemyPos.z,world->getParticles().back()));
                        break;
                    case 1:
                        enemies.push_back(new EnemyLVL3(enemyPos.x,3,enemyPos.z,world->getParticles().back()));
                        break;
                }
                
            }
        }
            
    }
    
    void spawnLVL7Enemies(int quantity) {
        if (enemies.size() == 0) {
            for (int i = 0; i < quantity; i++) {
                world->getParticles().push_back(new Particle());
                Vector3 enemyPos = getSafePosition();
                switch (i % 3) {
                    case 0:
                        enemies.push_back(new EnemyLVL1(enemyPos.x,3,enemyPos.z,world->getParticles().back()));
                        break;
                    case 1:
                        enemies.push_back(new EnemyLVL2(enemyPos.x,3,enemyPos.z,world->getParticles().back()));
                        break;
                    case 2:
                        enemies.push_back(new EnemyLVL3(enemyPos.x,3,enemyPos.z,world->getParticles().back()));
                        break;
                }
                
            }
        }
            
    }
    void spawnLVL8Enemies(int quantity) {
        if (enemies.size() == 0) {
            for (int i = 0; i < quantity; i++) {
                world->getParticles().push_back(new Particle());
                Vector3 enemyPos = getSafePosition();
                switch (i % 4) {
                    case 0:
                        enemies.push_back(new EnemyLVL1(enemyPos.x,3,enemyPos.z,world->getParticles().back()));
                        break;
                    case 1:
                        enemies.push_back(new EnemyLVL2(enemyPos.x,3,enemyPos.z,world->getParticles().back()));
                        break;
                    case 2:
                        enemies.push_back(new EnemyLVL3(enemyPos.x,3,enemyPos.z,world->getParticles().back()));
                        break;
                    case 3:
                        enemies.push_back(new EnemyLVL4(enemyPos.x,3,enemyPos.z,world->getParticles().back()));
                        break;
                }
                
            }
        }
    }
    
    
    /**
     * Odstrani vsetkych nepriatelov z hracej plochy
     **/
    void removeEnemies() {
        vector<Enemy*>::iterator en = enemies.begin();
        vector<Particle*> *worldParticles = &world->getParticles();
        for (; en != enemies.end(); en++) {
            vector<Particle*>::iterator index = find(worldParticles->begin(), worldParticles->end(), (*en)->particleRef);
            if (index != worldParticles->end()) {
                worldParticles->erase(index);
            }
        }
        enemies.clear();
    }
    
    void spawnPoints() {
        if (points.size() == 0) {
            for (int i = 0; i < gamePoints; i++) {
                points.push_back(new GamePoint());
            }
        }
    }
    
    void removePoints() {
        points.clear();
    }
    
    
public:
    /**
     * Inicializuje hraciu plochu a vygeneruje na nej 5 bodov (ktore sa neustale refreshuju ak niektory z nich zoberieme)
     **/
    Playground(GameCharacter *gmChr, Camera *cmr, ParticleWorld *world): gamechar(gmChr), camera(cmr), score(0), world(world) {

    }
    
    void draw() {
        drawFloor();
        gamechar->draw();
        drawPoints();
        drawEnemies();
        
    }
    void update() {
        // updat ukotvenia kamery
        
        if(playingGame){
            camera->cameraLink->changeAnchorPos(*gamechar->charPos);
        }
        
        updatePoints();
        updateEnemies();
        
        
        switch (score / levelChangeFactor) {
            case 1:
                spawnLVL1Enemies(10);
                break;
            case 2:
                spawnLVL2Enemies(20);
                break;
            case 3:
                spawnLVL3Enemies(20);
                break;
            case 4:
                spawnLVL4Enemies(10);
                break;
            case 5:
                spawnLVL5Enemies(15);
                break;
            case 6:
                spawnLVL6Enemies(15);
                break;
            case 7:
                spawnLVL7Enemies(18);
                break;
            case 8:
                spawnLVL8Enemies(20);
                break;
            case 9:
                goToMenu();
                break;
            default:
                if(!playingGame) {
                    // while in menu and not playing game
                    spawnLVL4Enemies(10);
                }
                break;
        }
        
        
        // update pozicie postavy
        gamechar->update();
        
        if (gamechar->lives <= 0) {
            goToMenu();
        }
    }
    
    /**
     * Obnovi hru na zaciatocne nastavenia - vsetko vynuluje
     **/
    void resetGame() {
        gamechar->respawn();
        camera->respawn();
        gameDuration = 0;
        score = 0;
        gamechar->lives = 3;
        endGame = false;
        removeEnemies();
        playingGame = true;
        inMenu = false;
        spawnPoints();
    }
    
    unsigned getEnemiesCount() {
        return (unsigned)enemies.size();
    }
    
    void superPower() {
        vector<Enemy*>::iterator en = enemies.begin();
        for (; en != enemies.end(); en++) {
            Vector3 jumpVel = (*en)->particleRef->getVelocity();
            jumpVel.addScaledVector(Vector3(0,1,0), 10.0f);
            (*en)->particleRef->setVelocity(jumpVel);
        }
    }
    
    void goToMenu() {
        if (!inMenu) {
            inMenu = true;
            endGame = true;
            playingGame = false;
            removeEnemies();
            spawnLVL4Enemies(15);
            gamechar->lives = 3;
            removePoints();
        }
    }
    
    void drawStats() {

        glDisable(GL_DEPTH_TEST);
        glMatrixMode( GL_PROJECTION ) ;
        glPushMatrix() ;
        glLoadIdentity();
        glMatrixMode( GL_MODELVIEW ) ;
        glPushMatrix() ;
        glLoadIdentity() ;
        
        
        
        if (!playingGame) {
            
            glColor4f(0.0, 0.0, 0.0, 0.8);
            glBegin(GL_QUADS);
            glVertex2f(-0.5, 0.7);
            glVertex2f(0.5, 0.7);
            glVertex2f(0.5, 0.9);
            glVertex2f(-0.5, 0.9);
            glEnd();
            
            glPushMatrix();
            glBegin(GL_QUADS);
            glVertex2f(-0.8, -0.6);
            glVertex2f(-0.4, -0.6);
            glVertex2f(-0.4, 0.6);
            glVertex2f(-0.8, 0.6);
            glEnd();
            glTranslatef(1.2, 0.0, 0.0);
            glBegin(GL_QUADS);
            glVertex2f(-0.8, -0.6);
            glVertex2f(-0.4, -0.6);
            glVertex2f(-0.4, 0.6);
            glVertex2f(-0.8, 0.6);
            glEnd();
            glPopMatrix();
            
            glBegin(GL_QUADS);
            glVertex2f(-0.2, -0.1);
            glVertex2f(0.2, -0.1);
            glVertex2f(0.2, 0.1);
            glVertex2f(-0.2, 0.1);
            glEnd();
            
        } else {
            glColor4f(1.0, 0.0, 0.0, 0.8);
            glBegin(GL_QUADS);
            
            glVertex2f(0.9, 0.8);
            glVertex2f(0.99, 0.8);
            glVertex2f(0.99, 0.99);
            glVertex2f(0.9, 0.99);
            glEnd();
            glColor4f(1.0, 0.5, 1.0, 0.8);
            for (int i = 0; i < gamechar->lives; i++) {
                glBegin(GL_QUADS);
                glVertex2f(-0.99+(float)i/10, 0.9);
                glVertex2f(-0.9+(float)i/10, 0.9);
                glVertex2f(-0.9+(float)i/10, 0.99);
                glVertex2f(-0.99+(float)i/10, 0.99);
                glEnd();
            }
        }
        
        glMatrixMode( GL_PROJECTION ) ;
        glPopMatrix() ;
        glMatrixMode( GL_MODELVIEW ) ;
        glPopMatrix() ;
        
        glColor3fv(WHITE);
        
        if (playingGame) {
            Draw::renderBitmapString(0.94, 0.88, GLUT_BITMAP_HELVETICA_18, convertInt(score));
        } else {
            Draw::renderBitmapString(-0.18, 0.78, GLUT_BITMAP_HELVETICA_18, "Bounce unlimited");
            Draw::renderBitmapString(-0.75, 0.5, GLUT_BITMAP_HELVETICA_18, "Controls:");
            Draw::renderBitmapString(-0.75, 0.4, GLUT_BITMAP_HELVETICA_18, "W-S-A-D = Movement");
            Draw::renderBitmapString(-0.75, 0.3, GLUT_BITMAP_HELVETICA_18, "Space = Jump");
            Draw::renderBitmapString(-0.75, 0.2, GLUT_BITMAP_HELVETICA_18, "P = super power!");
            Draw::renderBitmapString(-0.75, 0.1, GLUT_BITMAP_HELVETICA_18, "--------------");
            Draw::renderBitmapString(-0.75, 0.0, GLUT_BITMAP_HELVETICA_18, "R = start/reset game");
            Draw::renderBitmapString(-0.75, -0.1, GLUT_BITMAP_HELVETICA_18, "I = instructions");
            Draw::renderBitmapString(-0.75, -0.2, GLUT_BITMAP_HELVETICA_18, "ESC = quit game");

            Draw::renderBitmapString(0.45, 0.5, GLUT_BITMAP_HELVETICA_18, "About game:");
            Draw::renderBitmapString(0.45, 0.4, GLUT_BITMAP_HELVETICA_18, "You have 3 lives.");
            Draw::renderBitmapString(0.45, 0.3, GLUT_BITMAP_HELVETICA_18, "Dodge the enemies.");
            Draw::renderBitmapString(0.45, 0.2, GLUT_BITMAP_HELVETICA_18, "Reach the highest");
            Draw::renderBitmapString(0.45, 0.1, GLUT_BITMAP_HELVETICA_18, "score possible.");
            Draw::renderBitmapString(0.45, 0.0, GLUT_BITMAP_HELVETICA_18, "Every 5 points");
            Draw::renderBitmapString(0.45, -0.1, GLUT_BITMAP_HELVETICA_18, "enemies level up.");
            Draw::renderBitmapString(0.45, -0.2, GLUT_BITMAP_HELVETICA_18, "----------------");
            Draw::renderBitmapString(0.45, -0.3, GLUT_BITMAP_HELVETICA_18, "Autor: Michal Gallovic");
            
            std::stringstream scoreString;
            if (score/levelChangeFactor >8) {
                scoreString << "Congratulations! You have finished the game!";
                Draw::renderBitmapString(-0.155, -0.01, GLUT_BITMAP_HELVETICA_18, scoreString.str());
            } else {
                scoreString << "Your score: " << convertInt(score);
                if (score > 0) {
                    Draw::renderBitmapString(-0.12, -0.01, GLUT_BITMAP_HELVETICA_18, scoreString.str());
                } else {
                    Draw::renderBitmapString(-0.12, -0.01, GLUT_BITMAP_HELVETICA_18, "PRESS R TO PLAY");
                }
            }
            
            
        }
        
        
        glEnable(GL_DEPTH_TEST);
        

    }
    
    
};

class PhysicsSimulation : public MassAggregateApplication {
protected:
    map<int, real> keysPressed;
    GameCharacter *gamechar;
    Camera *camera;
    WallContacts wallContactGenerator;
    EnemyContacts enemyContactGenerator;
    Playground *playground;
public:
    PhysicsSimulation();
    virtual ~PhysicsSimulation();
    virtual void update();
    virtual void display();
    virtual void initGraphics();
    virtual const char *getTitle();
    virtual void key(int key);
    virtual void keyUP(int key);
    
};





PhysicsSimulation::PhysicsSimulation():MassAggregateApplication(DEFAULT_PARTICLES) {
    keysPressed = {{'w',0.0f},{'s',0.0f},{'a',0.0f},{'d',0.0f}};
    gamechar = new GameCharacter(&particleArray[GAME_CHAR_INDEX]);
    camera = new Camera(&particleArray[CAMERA_INDEX]);
    playground = new Playground(gamechar, camera, &world);
    playground->gameDuration = 0;
    for (int i = 0; i < DEFAULT_PARTICLES; i++) {
        particleArray[i].clearAccumulator();
    }
    
    wallContactGenerator.init(&world.getParticles());
    world.getContactGenerators().push_back(&wallContactGenerator);
    enemyContactGenerator.init(&world.getParticles());
    world.getContactGenerators().push_back(&enemyContactGenerator);
    //camera<---|gamechar link
    camera->cameraLink = new CameraSpring(new Vector3(25,25,25), 3.0, 1.0);
    world.getForceRegistry().add(camera->particleRef, camera->cameraLink);
    
    
}

void PhysicsSimulation::initGraphics() {
    MassAggregateApplication::initGraphics();
    
    
    //osvetlenie
    //ambient
    GLfloat ambientColor[] = {0.5f, 0.5f, 0.5f, 1.0f};
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambientColor);

    //positioned
    GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat mat_shininess[] = { 50.0 };
    GLfloat mat_emission[] = {0.1,0.1,0.1,1.0};
    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
    glMaterialfv(GL_FRONT, GL_EMISSION, mat_emission);
    
    GLfloat lightColor0[] = {0.5f, 0.5f, 0.5f,1.0f};
    GLfloat lightPos0[] = {0.0f,1.0f,0.0f,0.0f};
    GLfloat lightAmbient0[] = { 0.35, 0.35, 0.35, 1.0 };
    glLightfv(GL_LIGHT0, GL_DIFFUSE, lightColor0);
    glLightfv(GL_LIGHT0, GL_POSITION, lightPos0);
    glLightfv(GL_LIGHT0, GL_AMBIENT, lightAmbient0);

    
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable (GL_BLEND);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    glEnable(GL_NORMALIZE);
}

void PhysicsSimulation::update() {
    if (playground->endGame) {
        camera->cameraLink->changeAnchorPos(Vector3(25,25,25));
    }
        // update sceny
        playground->update();
        MassAggregateApplication::update();
        playground->gameDuration += TimingData::get().lastFrameDuration;
        
        // pohyb postavy
        gamechar->move(keysPressed, camera);
        
        world.setMaxContacts(playground->getEnemiesCount()*5 + DEFAULT_PARTICLES);
}



void PhysicsSimulation::display() {
    // Clear the view port and set the camera direction
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();

    Vector3 pos = gamechar->getPosition();
    Vector3 cPos = camera->particleRef->getPosition();
    
    
    
    gluLookAt(cPos.x, cPos.y, cPos.z,  pos.x, pos.y, pos.z,  0.0, 1.0, 0.0);
    playground->draw();
    playground->drawStats();
    

}

void PhysicsSimulation::key(int key) {
    // zvacsi smer pohybu o delta pre novo stlacenu klavesu
    real delta = 0.1f;
    if (playingGame) {
        switch(key)
        {
            case 'w':
                keysPressed['w'] += delta;
                break;
            case 's':
                keysPressed['s'] += delta;
                break;
            case 'a':
                keysPressed['a'] += delta;
                break;
            case 'd':
                keysPressed['d'] += delta;
                break;
            case 'p':
                playground->superPower();
                break;
            case SPACE_KEY:
                gamechar->jump();
                break;
            default:
                MassAggregateApplication::key(key);
        }
    }
    
    switch (key) {
        case 'r':
            // reset game
            playground->resetGame();
            break;
        case 'i':
            // do menu
            playground->goToMenu();
            break;
        case 27:
            exit(1);
            break;
    }
    
    
    map<int,real>::iterator i = keysPressed.begin();
    // + pre diagonalne pohyby - zvacsi o delta aj klavesu, ktoru sme nepustili
    for (; i != keysPressed.end(); i++) {
        if (i->second != 0.0f && i->first != key) {
            i->second += delta;
        }
    }
    
    
    
}

void PhysicsSimulation::keyUP(int key) {
    // smer pohybu, ktory sme prestali zvacsovat automaticky nastav na 0
    // zrychlenie sa v tomto smere nastavi na 0, ale rychlost ostane
    // a bude sa postupne zmensovat
    switch(key)
    {
        case 'w':
            keysPressed['w'] = 0.0f;
            break;
        case 's':
            keysPressed['s'] = 0.0f;
            break;
        case 'a':
            keysPressed['a'] = 0.0f;
            break;
        case 'd':
            keysPressed['d'] = 0.0f;
            break;
    }
}

const char* PhysicsSimulation::getTitle() {
    return "Bounce unlimited - Michal Gallovic";
}

PhysicsSimulation::~PhysicsSimulation() {
    
}

Application* getApplication() {
    return new PhysicsSimulation();
}