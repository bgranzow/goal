#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apf.h>
#include <PCU.h>
#include <apfShape.h>
#include <mthQR.h>
#include <set>
#include <pcu_util.h> // for assert
#include <iostream>

namespace enrich {
/* overall information useful during recovery */
struct Recovery {
  apf::Mesh* mesh;
  /* mesh dimension, so far handling 2 and 3 */
  int dim;
  /* order of output field and fit polynomial.
	 determined by the order of the mesh's
     coordinate field. */
  int order;
  /* the number of terms in a polynomial of the above order
     defined over ND space: f(x,y,z) in 3D */
  int polynomial_terms;
  /* the number of integration points per element.
     this immediately assumes the mesh has one element type.
     warning: the number of points per element is not necessarily
     determined by the order of the output field !
     users sometimes give us "too much" information, i.e. 5 integration
     points per tet to derive a linear field. this is ok, the
     way this is programmed it should handle any nonzero number of points
     per element regardless of output order (by growing bigger patches) */
  int points_per_element;
  /* input field containing integration point data for all elements */
  apf::Field* f;
  /* output field containing enriched nodal data */
  apf::Field* f_star;
};

static int determinePointsPerElement(apf::Field* f)
{
  apf::Mesh* m = apf::getMesh(f);
  int element_type = apf::getFirstType(m, m->getDimension());
  apf::FieldShape* s = apf::getShape(f);
  return s->countNodesOn(element_type);
}

static apf::Field* makeRecoveredField(Recovery* r)
{
  std::string name = "enrich_";
  name += apf::getName(r->f);
  return apf::createLagrangeField(
        r->mesh, name.c_str(), apf::getValueType(r->f), r->order);
}

static int countPolynomialTerms(int dim, int order)
{
  switch (dim) {
    case 2:
      return ((order + 1) * (order + 2)) / 2;
    case 3:
      return ((order + 1) * (order + 2) * (order + 3)) / 6;
    default:
      apf::fail("bad dim in countPolynomialTerms");
      return -1;
  }
}

static void setupRecovery(Recovery* r, apf::Field* f)
{
  r->mesh = apf::getMesh(f);
  r->dim = r->mesh->getDimension();
  r->order = r->mesh->getShape()->getOrder();
  r->polynomial_terms = countPolynomialTerms(r->dim, r->order);
  r->points_per_element = determinePointsPerElement(f);
  r->f = f;
  r->f_star = makeRecoveredField(r);
}

struct Samples {
  Samples():num_points(0) {}
  void allocate(int np, int nc)
  {
    num_points = np;
    points.allocate(np);
    values.allocate(np);
    for (int i=0; i < np; ++i)
      values[i].allocate(nc);
  }
  int num_points;
  apf::NewArray<apf::Vector3> points;
  apf::NewArray<apf::NewArray<double> > values;
};

struct QRDecomp {
  mth::Matrix<double> Q;
  mth::Matrix<double> R;
};

typedef std::set<apf::MeshEntity*> EntitySet;

struct Patch {
  apf::Mesh* mesh;
  Recovery* recovery;
  /* the entity around which the patch
     is centered. a patch collects elements
     around this entity and then their integration
     points will be used to recover values for
     all nodes on this entity */
  apf::MeshEntity* entity;
  EntitySet elements;
  Samples samples;
  QRDecomp qr;
};

static void setupPatch(Patch* p, Recovery* r)
{
  p->mesh = r->mesh;
  p->recovery = r;
  p->entity = 0;
}

static void startPatch(Patch* p, apf::MeshEntity* e)
{
  p->elements.clear();
  p->entity = e;
}

static void addElementToPatch(Patch* p, apf::MeshEntity* e)
{
  p->elements.insert(e);
}

static void addElementsToPatch(Patch* p, apf::DynamicArray<apf::MeshEntity*>& es)
{
  for (std::size_t i=0; i < es.getSize(); ++i)
	{
    		addElementToPatch(p, es[i]);
	}
}



static bool getInitialPatch(Patch* p)
{
  apf::DynamicArray<apf::MeshEntity*> adjacent;
  p->mesh->getAdjacent(p->entity, p->recovery->dim, adjacent);
  addElementsToPatch(p, adjacent);
  return true;
}



std::vector<apf::Vector3> createPoints(int dim, int n)
{
	// do the local coords range from 0-1?
	std::vector<apf::Vector3> pointCoords;
	if (dim == 2)
	{
		pointCoords.push_back(apf::Vector3(.333,.333,0));
		pointCoords.push_back(apf::Vector3(.1,.1,0));
		pointCoords.push_back(apf::Vector3(.9,.1,0));
		pointCoords.push_back(apf::Vector3(.1,.9,0));
		pointCoords.push_back(apf::Vector3(.05,.7,0));
		pointCoords.push_back(apf::Vector3(.25,.7,0));
		pointCoords.push_back(apf::Vector3(.15,.3,0));
		pointCoords.push_back(apf::Vector3(.7,.05,0));
		pointCoords.push_back(apf::Vector3(.7,.25,0));
		pointCoords.push_back(apf::Vector3(.3,.15,0));
	}
	if (dim == 3)
	{
		pointCoords.push_back(apf::Vector3(.05,.05,.05));
		pointCoords.push_back(apf::Vector3(.05,.9,.05));
		pointCoords.push_back(apf::Vector3(.05,.05,.9));
		pointCoords.push_back(apf::Vector3(.9,.05,.05));
		pointCoords.push_back(apf::Vector3(0,.33,.33));
		pointCoords.push_back(apf::Vector3(.33,0,.33));
		pointCoords.push_back(apf::Vector3(.33,.33,0));
		pointCoords.push_back(apf::Vector3(.33,.33,.33));
		pointCoords.push_back(apf::Vector3(.1,.05,.05));
		pointCoords.push_back(apf::Vector3(.05,.1,.05));
		pointCoords.push_back(apf::Vector3(.05,.05,.1));
		pointCoords.push_back(apf::Vector3(.4,.4,.1));
		pointCoords.push_back(apf::Vector3(.4,.1,.4));
		pointCoords.push_back(apf::Vector3(.1,.4,.4));
		pointCoords.push_back(apf::Vector3(.6,.1,.1));
		pointCoords.push_back(apf::Vector3(.1,.6,.1));
		pointCoords.push_back(apf::Vector3(.1,.1,.6));
		pointCoords.push_back(apf::Vector3(.15,.25,.35));
		pointCoords.push_back(apf::Vector3(.35,.15,.25));
		pointCoords.push_back(apf::Vector3(.25,.35,.15));
	}
	pointCoords.resize(n);
	return pointCoords;
}



static void getSamplePoints(Patch* p)
{
	Recovery* r = p->recovery;
	Samples* s = &p->samples;
	int np = countPolynomialTerms(r->dim, r->order);

	// use all the samples	
	np=countPolynomialTerms(r->dim, 3)*(p->elements.size());

	int nc = apf::countComponents(r->f);
	s->allocate(np,nc); // does both points and values
	int samplesPerElement=int(ceil(float(np)/float(p->elements.size())));

	// use all the samples:
	samplesPerElement=countPolynomialTerms(r->dim, r->order); // only works if you want cubics (need more neutral points for quad)

	//std::cout<<"dim: "<<(r->dim)<<" order: "<<(r->order)<<" np: "<<np<<" pp: "<<(p->elements.size())<<" samples per element: "<<samplesPerElement<<std::endl;
  
	std::vector<apf::Vector3> pointCoords = createPoints(r->dim, samplesPerElement);
	std::size_t i = 0;
	APF_ITERATE(EntitySet, p->elements, it) 
	{
		apf::MeshElement* me = apf::createMeshElement(r->mesh, *it);
		apf::Element* el=apf::createElement(r->f,me);
		for (int l = 0; l < samplesPerElement; ++l) 
		{
			if (i>=size_t(np))
				break;
			//std::cout<<i<<" "<<size_t(np-1)<<" "<<l<<std::endl;
			apf::getComponents(el, pointCoords[l], &(s->values[i][0]));
      			apf::mapLocalToGlobal(me, pointCoords[l], s->points[i]);
			++i;
      			
    		}
    		apf::destroyMeshElement(me);
		apf::destroyElement(el);
  	}
}


static void evalPolynomialTerms(
    int dim, int order,
    apf::Vector3 const& point,
    mth::Vector<double>& terms)
{
  apf::Vector3 const& x = point;
  switch (dim) {
  case 2:
    switch (order) {
    case 1:
      terms.resize(3);
      terms(0) = 1.0;
      terms(1) = x[0];
      terms(2) = x[1];
      return;
    case 2:
      terms.resize(6);
      terms(0) = 1.0;
      terms(1) = x[0];
      terms(2) = x[1];
      terms(3) = x[0]*x[1];
      terms(4) = x[0]*x[0];
      terms(5) = x[1]*x[1];
      return;
	case 3: //added 2D cubic
      terms.resize(10);
      terms(0) = 1.0;
      terms(1) = x[0];
      terms(2) = x[1];
      terms(3) = x[0]*x[1];
      terms(4) = x[0]*x[0];
      terms(5) = x[1]*x[1];
      terms(6) = x[0]*x[0]*x[0];
      terms(7) = x[0]*x[0]*x[1];
      terms(8) = x[1]*x[1]*x[0];
      terms(9) = x[1]*x[1]*x[1];
      return;
    default:
      apf::fail("ENRICH: invalid 2D polynomial order");
    }
  case 3:
    switch (order) {
    case 1:
      terms.resize(4);
      terms(0) = 1.0;
      terms(1) = x[0];
      terms(2) = x[1];
      terms(3) = x[2];
      return;
    case 2:
      terms.resize(10);
      terms(0) = 1.0;
      terms(1) = x[0];
      terms(2) = x[1];
      terms(3) = x[2];
      terms(4) = x[0]*x[1];
      terms(5) = x[1]*x[2];
      terms(6) = x[2]*x[0];
      terms(7) = x[0]*x[0];
      terms(8) = x[1]*x[1];
      terms(9) = x[2]*x[2];
      return;
	case 3: // added 3D cubic
      terms.resize(20);
      terms(0) = 1.0;
      terms(1) = x[0];
      terms(2) = x[1];
      terms(3) = x[2];
      terms(4) = x[0]*x[1];
      terms(5) = x[1]*x[2];
      terms(6) = x[2]*x[0];
      terms(7) = x[0]*x[0];
      terms(8) = x[1]*x[1];
      terms(9) = x[2]*x[2];
	  terms(10) = x[0]*x[0]*x[0];
      terms(11) = x[0]*x[0]*x[1];
      terms(12) = x[0]*x[0]*x[2];
      terms(13) = x[1]*x[1]*x[0];
      terms(14) = x[1]*x[1]*x[1];
      terms(15) = x[1]*x[1]*x[2];
	  terms(16) = x[2]*x[2]*x[0];
      terms(17) = x[2]*x[2]*x[1];
      terms(18) = x[2]*x[2]*x[2];
      terms(19) = x[0]*x[1]*x[2];
      return;
    default:
      apf::fail("ENRICH: invalid 3D polynomial order");
    }
  default:
    apf::fail("ENRICH: invalid polynomial order");
  }
}


static bool preparePolynomialFit(
    int dim,
    int order,
    int num_points,
    apf::NewArray<apf::Vector3> const& points,
    QRDecomp& qr)
{
  unsigned m = num_points;
  unsigned n = countPolynomialTerms(dim, order);
  PCU_ALWAYS_ASSERT(m >= n);
  mth::Matrix<double> A(m,n);
  mth::Vector<double> p;
  for (unsigned i = 0; i < m; ++i) {
    evalPolynomialTerms(dim, order, points[i], p);
    for (unsigned j = 0; j < p.size(); ++j)
      A(i,j) = p(j);
  }
///////////std::cout<<A<<std::endl;
  unsigned rank = mth::decomposeQR(A, qr.Q, qr.R);
  return rank == A.cols();
}


static void runPolynomialFit(QRDecomp const& qr,
                             mth::Vector<double> const& values,
                             mth::Vector<double>& coeffs)
{
  mth::solveFromQR(qr.Q, qr.R, values, coeffs);
}


static double evalPolynomial(int dim, int order, apf::Vector3& point,
    mth::Vector<double>& coeffs)
{
  mth::Vector<double> terms;
  evalPolynomialTerms(dim, order, point, terms);
  return coeffs * terms;
}


static void runEnrich(Patch* p)
{
  Recovery* r = p->recovery;
  apf::Mesh* m = r->mesh;
  Samples* s = &p->samples;
  //getSampleValues(p); // NOT ANYMORE HAH
  int num_components = apf::countComponents(r->f_star);
  int num_nodes = m->getShape()->countNodesOn(m->getType(p->entity));
  mth::Vector<double> values(s->num_points);
  apf::NewArray<apf::Vector3> nodal_points(num_nodes);
  apf::NewArray<apf::NewArray<double> > recovered_values(num_nodes);
  for (int i = 0; i < num_nodes; ++i) {
    recovered_values[i].allocate(num_components);
    m->getPoint(p->entity, i, nodal_points[i]);
  }
  for (int i = 0; i < num_components; ++i) {
    for (int j = 0; j < s->num_points; ++j)
      values(j) = s->values[j][i];
    mth::Vector<double> coeffs;
    runPolynomialFit(p->qr, values, coeffs);
    for (int j = 0; j < num_nodes; ++j)
      recovered_values[j][i] = evalPolynomial(
          r->dim, r->order, nodal_points[j], coeffs);
  }
  for (int i = 0; i < num_nodes; ++i)
    apf::setComponents(r->f_star, p->entity, i, &(recovered_values[i][0]));
}


static bool buildPatch(Patch* p)
{
  if (!getInitialPatch(p)) return false;
  return true;
}


apf::Field* enrichField(apf::Field* f, int desiredOrder)
{
	Recovery recovery; // old stuff
	setupRecovery(&recovery, f);

	recovery.order=desiredOrder; 

	Patch patch;
	apf::Mesh* m = apf::getMesh(f); // get mesh
	for (int d=0; d <= (m->getDimension()-1); ++d) // loop over dims
	{
		apf::MeshEntity* s;
		apf::MeshIterator* it = m->begin(d); // everything but highest dim
		while ((s = m->iterate(it))) // loop over entities of that dim
		{
			setupPatch(&patch, &recovery); // create patch of elements adjacent to this mesh entity
			startPatch(&patch, s);
			buildPatch(&patch);
			getSamplePoints(&patch); // some sample points, int points and some extras
			preparePolynomialFit(recovery.dim, recovery.order, patch.samples.num_points, patch.samples.points, patch.qr);
			runEnrich(&patch); // this will run the least squares and add result to new field
		}
		m->end(it);
	}
  return recovery.f_star;

}

}

