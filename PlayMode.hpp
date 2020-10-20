#include "Mode.hpp"

#include "Scene.hpp"
#include "WalkMesh.hpp"
#include "Sound.hpp"

#include <glm/glm.hpp>

#include <vector>
#include <deque>

struct PlayMode : Mode {
	PlayMode();
	virtual ~PlayMode();

	//functions called by main loop:
	virtual bool handle_event(SDL_Event const &, glm::uvec2 const &window_size) override;
	virtual void update(float elapsed) override;
	virtual void draw(glm::uvec2 const &drawable_size) override;

	//----- game state -----

	//input tracking:
	struct Button {
		uint8_t downs = 0;
		uint8_t pressed = 0;
	} left, right, down, up;

	//local copy of the game scene (so code can change it during gameplay):
	Scene scene;

	std::vector<std::string> hints;

	// infinity asset movement
	Scene::Transform *infinity = nullptr;
	float total_inf_translation = 0.0f;
	float inf_speed = 0.1f;

	// lines asset movement
	Scene::Transform *line1 = nullptr;
	Scene::Transform *line2 = nullptr;
	Scene::Transform *line3 = nullptr;
	glm::quat line1_base_rotation;
	glm::quat line2_base_rotation;
	glm::quat line3_base_rotation;
	float swing = 0.0f;

	Scene::Transform *box = nullptr;
	float box_speed = 1.0f;
	float cone_speed = -2.0f;
	float up_and_down = 0.04f;
	float cone_total = 2.0f;

	Scene::Transform *cone1 = nullptr;
	Scene::Transform *cone2 = nullptr;
	Scene::Transform *cylinder_long = nullptr;

	float cylinder_long_total = 0.0f;
	
	Scene::Transform *ball1 = nullptr;
	Scene::Transform *ball2 = nullptr;
	bool try_to_catch = false;
	bool ball1_caught = false;
	bool ball2_caught = false;
	bool reset = false;

	glm::vec3 turning_point = glm::vec3(1.59654f, -0.147706f, 6.04984f);

	int cur_string_index = 0;
	bool to_do = true;
	float update_window = 1.0f;

	//player info:
	struct Player {
		WalkPoint at;
		//transform is at player's feet and will be yawed by mouse left/right motion:
		Scene::Transform *transform = nullptr;
		//camera is at player's head and will be pitched by mouse up/down motion:
		Scene::Camera *camera = nullptr;
	} player;
};
