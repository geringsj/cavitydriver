#include "CameraSystem.h"

/*! 
 * constructor creates a camera system with the camera in position (0,0,0)
 * and up-vector (0,1,0) forward-vector (0,0,1) and right-vector (1,0,0).
 */
CameraSystem::CameraSystem()
{
	cam_pos = glm::vec3(0.0f,0.0f,0.0f);
	up_vector = glm::vec3(0.0f,1.0f,0.0f);
	forward_vector = glm::vec3(0.0f,0.0f,1.0f);
	right_vector = glm::vec3(1.0f,0.0f,0.0f);
	center = glm::vec3(0.0f,0.0f,0.0f);
	translation = glm::mat4(1.0f);
	rotation = glm::mat4(1.0f);
}

/*! 
 * constructor creates a camera system with the given parameters.
 */
CameraSystem::CameraSystem(glm::vec3 p_cam_pos,glm::vec3 p_up_vector,glm::vec3 p_forward_vector,glm::vec3 p_right_vector)
{
	cam_pos = p_cam_pos;
	up_vector = p_up_vector;
	forward_vector = p_forward_vector;
	right_vector = p_right_vector;
	center = glm::vec3(0.0f,0.0f,0.0f);
	translation = glm::mat4(1.0f);
	rotation = glm::mat4(1.0f);
}

/*!
 * Calculate the rotation with the given vector and the angle alpha.
 * The rotation affects all vectors (up,right and forward). Also the
 * center is affected.
 */
void CameraSystem::Rotation(glm::vec3 v,float alpha)
{
	rotation = glm::rotate(glm::mat4(1.0f),alpha,v);
	up_vector = MakeVec4ToVec3(rotation * MakeVec3ToVec4(up_vector));
	forward_vector = MakeVec4ToVec3(rotation * MakeVec3ToVec4(forward_vector));
	right_vector = MakeVec4ToVec3(rotation * MakeVec3ToVec4(right_vector));
	center = MakeVec4ToVec3(rotation * MakeVec3ToVec4(center));
}

/*!
 * Calculate the translation in direction of vector t with the step
 * size. The camera position and the center vector are affected.
 */
void CameraSystem::Translation(glm::vec3 t, float step_size)
{
	translation = glm::translate(glm::mat4(1.0f),step_size*t);
	cam_pos = MakeVec4ToVec3(translation * MakeVec3ToVec4(cam_pos));
	center = MakeVec4ToVec3(translation * MakeVec3ToVec4(center));
}

/*!
 * Create the view matrix by using glm lookAt. Use the position, center and
 * up-vector. Return the 4x4 matrix.
 */
glm::mat4 CameraSystem::GetViewMatrix()
{
	return glm::lookAt(cam_pos,center,up_vector);
}

/*!
 * Convert a 4D vector to a 3D vector by dividing with the w component.
 */
glm::vec3 CameraSystem::MakeVec4ToVec3(glm::vec4 to_shrink)
{
	glm::vec3 tmp_vec;
	tmp_vec.x = to_shrink.x/to_shrink.w;
	tmp_vec.y = to_shrink.y/to_shrink.w;
	tmp_vec.z = to_shrink.z/to_shrink.w;
	return tmp_vec;
}

/*!
 * Convert a 3D vector to a 4D vector by adding a w component with value 1.
 */
glm::vec4 CameraSystem::MakeVec3ToVec4(glm::vec3 to_grow)
{
	return glm::vec4(to_grow,1.0f);
}

CameraSystem::~CameraSystem()
{

}