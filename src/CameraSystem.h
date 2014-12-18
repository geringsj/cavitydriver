#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
//! GLM includes
#include "glm/glm.hpp"
#include "glm/gtc/type_ptr.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtx/rotate_vector.hpp"

class CameraSystem
{
public:
	//! constructor creates a camera system with the camera in position (0,0,0)
	//! and up-vector (0,1,0) forward-vector (0,0,1) and right-vector (1,0,0)
	CameraSystem();
	//! constructor creates a camera system with the given parameters
	CameraSystem(glm::vec3 p_cam_pos,glm::vec3 p_up_vector,glm::vec3 p_forward_vector,glm::vec3 p_right_vector);
	//! rotation around (up,forward,right)-vector v with angle alpha (in degree)
	void Rotation(glm::vec3 v,float alpha);
	//! translation with a step in (up,forward,right)-vector direction
	void Translation(glm::vec3 t, float step_size);
	//! return the view matrix
	glm::mat4 GetViewMatrix();
	//! returns the camera position
	glm::vec3 GetCamPos() { return cam_pos; }
	//! returns the up-vector
	glm::vec3 GetUpVector() { return up_vector; }
	//! returns the forwards vector
	glm::vec3 GetForwardVector() { return forward_vector; }
	//! returns the right vector
	glm::vec3 GetRightVector() { return right_vector; }
	//! returns the center vector
	glm::vec3 GetCenterVector() { return center; }
	~CameraSystem();
private:
	//! convert from 3D vector to a 4D vector by adding the w component with w = 1
	glm::vec3 MakeVec4ToVec3(glm::vec4 to_shrink);
	//! convert from 4D vector to a 3D vector by dividing by the w component
	glm::vec4 MakeVec3ToVec4(glm::vec3 to_grow);
	//! stores the camera position
	glm::vec3 cam_pos;
	//! stores the up-vector
	glm::vec3 up_vector;
	//! stores the forward vector
	glm::vec3 forward_vector;
	//! stores the right vector
	glm::vec3 right_vector;
	//! stores the center vector, i.e. the positin we look at in world space
	glm::vec3 center;
	//! stores the translation
	glm::mat4 translation;
	//! stores the rotation
	glm::mat4 rotation;
};