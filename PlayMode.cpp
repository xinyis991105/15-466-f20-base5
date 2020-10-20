#include "PlayMode.hpp"

#include "LitColorTextureProgram.hpp"

#include "DrawLines.hpp"
#include "Mesh.hpp"
#include "Load.hpp"
#include "gl_errors.hpp"
#include "data_path.hpp"

#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/quaternion.hpp>

#include <random>

GLuint phonebank_meshes_for_lit_color_texture_program = 0;
Load< MeshBuffer > phonebank_meshes(LoadTagDefault, []() -> MeshBuffer const * {
	MeshBuffer const *ret = new MeshBuffer(data_path("game5.pnct"));
	phonebank_meshes_for_lit_color_texture_program = ret->make_vao_for_program(lit_color_texture_program->program);
	return ret;
});

Load< Scene > phonebank_scene(LoadTagDefault, []() -> Scene const * {
	return new Scene(data_path("game5.scene"), [&](Scene &scene, Scene::Transform *transform, std::string const &mesh_name){
		Mesh const &mesh = phonebank_meshes->lookup(mesh_name);

		scene.drawables.emplace_back(transform);
		Scene::Drawable &drawable = scene.drawables.back();

		drawable.pipeline = lit_color_texture_program_pipeline;

		drawable.pipeline.vao = phonebank_meshes_for_lit_color_texture_program;
		drawable.pipeline.type = mesh.type;
		drawable.pipeline.start = mesh.start;
		drawable.pipeline.count = mesh.count;

	});
});

WalkMesh const *walkmesh = nullptr;
Load< WalkMeshes > phonebank_walkmeshes(LoadTagDefault, []() -> WalkMeshes const * {
	WalkMeshes *ret = new WalkMeshes(data_path("game5.w"));
	walkmesh = &ret->lookup("Plane");
	return ret;
});

Load< Sound::Sample > bgm_load(LoadTagDefault, []() -> Sound::Sample const * {
	return new Sound::Sample(data_path("bgm.wav"));
});

Load< Sound::Sample > sound_effect_load(LoadTagDefault, []() -> Sound::Sample const * {
	return new Sound::Sample(data_path("there_you_go.wav"));
});

PlayMode::PlayMode() : scene(*phonebank_scene) {
	//create a player transform:
	for (auto &transform : scene.transforms) {
		if (transform.name == "Infinity") {
			infinity = &transform;
		} else if (transform.name == "Lines") {
			line1 = &transform;
			line1_base_rotation = line1->rotation;
		} else if (transform.name == "Lines.001") {
			line2 = &transform;
			line2_base_rotation = line2->rotation;
		} else if (transform.name == "Lines.002") {
			line3 = &transform;
			line3_base_rotation = line3->rotation;
		} else if (transform.name == "Cube") {
			box = &transform;
		} else if (transform.name == "Cube.001") {
			ball1 = &transform;
		} else if (transform.name == "Cube.002") {
			ball2 = &transform;
		} else if (transform.name == "Cylinder") {
			cone1 = &transform;
		} else if (transform.name == "Cylinder.002") {
			cylinder_long = &transform;
		} else if (transform.name == "Cylinder.003") {
			cone2 = &transform;
		}
	}
	scene.transforms.emplace_back();
	player.transform = &scene.transforms.back();

	//create a player camera attached to a child of the player transform:
	scene.transforms.emplace_back();
	scene.cameras.emplace_back(&scene.transforms.back());
	player.camera = &scene.cameras.back();
	player.camera->fovy = glm::radians(60.0f);
	player.camera->near = 0.01f;
	player.camera->transform->parent = player.transform;

	//player's eyes are 1.8 units above the ground:
	player.camera->transform->position = glm::vec3(0.0f, 0.0f, 1.8f);

	//rotate camera facing direction (-z) to player facing direction (+y):
	player.camera->transform->rotation = glm::angleAxis(glm::radians(90.0f), glm::vec3(1.0f, 0.0f, 0.0f));

	//start player walking at nearest walk point:
	player.at = walkmesh->nearest_walk_point(player.transform->position);

	Sound::loop_3D(*bgm_load, 0.15f, player.transform->position, 10.0f);

	hints.push_back(std::string("Left click to control your angle."));
	hints.push_back(std::string("Press Escape to ungrab your mouse."));
	hints.push_back(std::string("Press W to move forward."));
	hints.push_back(std::string("Catch the two rotating balls by pressing C."));
	hints.push_back(std::string("Press R to restart."));
}

PlayMode::~PlayMode() {
}

bool PlayMode::handle_event(SDL_Event const &evt, glm::uvec2 const &window_size) {

	if (evt.type == SDL_KEYDOWN) {
		if (evt.key.keysym.sym == SDLK_ESCAPE) {
			if (cur_string_index == 1) to_do = false;
			SDL_SetRelativeMouseMode(SDL_FALSE);
			return true;
		} else if (evt.key.keysym.sym == SDLK_w) {
			if (cur_string_index == 2) to_do = false;
			up.downs += 1;
			up.pressed = true;
			return true;
		} else if (evt.key.keysym.sym == SDLK_c && cur_string_index == 3) {
			try_to_catch = true;
			return true;
		} else if (evt.key.keysym.sym == SDLK_r) {
			if (ball1_caught && ball2_caught) {
				reset = true;
				return true;
			}
		}
	} else if (evt.type == SDL_KEYUP) {
		if (evt.key.keysym.sym == SDLK_w) {
			up.pressed = false;
			return true;
		}
	} else if (evt.type == SDL_MOUSEBUTTONDOWN) {
		if (SDL_GetRelativeMouseMode() == SDL_FALSE) {
			if (cur_string_index == 0) to_do = false;
			SDL_SetRelativeMouseMode(SDL_TRUE);
			return true;
		}
	} else if (evt.type == SDL_MOUSEMOTION) {
		if (SDL_GetRelativeMouseMode() == SDL_TRUE) {
			glm::vec2 motion = glm::vec2(
				evt.motion.xrel / float(window_size.y),
				-evt.motion.yrel / float(window_size.y)
			);
			glm::vec3 up = walkmesh->to_world_smooth_normal(player.at);
			player.transform->rotation = glm::angleAxis(-motion.x * player.camera->fovy, up) * player.transform->rotation;

			float pitch = glm::pitch(player.camera->transform->rotation);
			pitch += motion.y * player.camera->fovy;
			//camera looks down -z (basically at the player's feet) when pitch is at zero.
			pitch = std::min(pitch, 0.95f * 3.1415926f);
			pitch = std::max(pitch, 0.05f * 3.1415926f);
			player.camera->transform->rotation = glm::angleAxis(pitch, glm::vec3(1.0f, 0.0f, 0.0f));

			return true;
		}
	}

	return false;
}

void PlayMode::update(float elapsed) {
	if (ball1_caught && ball2_caught) {
		if (reset) {
			player.camera->fovy = glm::radians(60.0f);
			player.camera->near = 0.01f;
			player.camera->transform->position = glm::vec3(0.0f, 0.0f, 1.8f);
			player.camera->transform->rotation = glm::angleAxis(glm::radians(90.0f), glm::vec3(1.0f, 0.0f, 0.0f));
			player.transform->position = glm::vec3(0.0f, 0.0f, 0.0f);
			player.at = walkmesh->nearest_walk_point(player.transform->position);
			ball1->scale = glm::vec3(1.0f);
			ball2->scale = glm::vec3(1.0f);
			total_inf_translation = 0.0f;
			inf_speed = 0.1f;
			swing = 0.0f;
			box_speed = 1.0f;
			cone_speed = -2.0f;
			up_and_down = 0.04f;
			cone_total = 2.0f;
			cylinder_long_total = 0.0f;
			try_to_catch = false;
			ball1_caught = false;
			ball2_caught = false;
			reset = false;
			cur_string_index = 3;
			to_do = true;
			update_window = 1.0f;
		} else {
			to_do = true;
			cur_string_index = 4;
		}
	} else {
		//player walking:
		{
			//combine inputs into a move:
			glm::vec3 old_pos = player.transform->position;
			WalkPoint old_at = player.at;
			constexpr float PlayerSpeed = 3.0f;
			glm::vec2 move = glm::vec2(0.0f);
			if (left.pressed && !right.pressed) move.x =-1.0f;
			if (!left.pressed && right.pressed) move.x = 1.0f;
			if (down.pressed && !up.pressed) move.y =-1.0f;
			if (!down.pressed && up.pressed) move.y = 1.0f;

			//make it so that moving diagonally doesn't go faster:
			if (move != glm::vec2(0.0f)) move = glm::normalize(move) * PlayerSpeed * elapsed;

			//get move in world coordinate system:
			glm::vec3 remain = player.transform->make_local_to_world() * glm::vec4(move.x, move.y, 0.0f, 0.0f);

			//using a for() instead of a while() here so that if walkpoint gets stuck in
			// some awkward case, code will not infinite loop:
			for (uint32_t iter = 0; iter < 10; ++iter) {
				if (remain == glm::vec3(0.0f)) break;
				WalkPoint end;
				float time;
				walkmesh->walk_in_triangle(player.at, remain, &end, &time);
				player.at = end;
				if (time == 1.0f) {
					//finished within triangle:
					remain = glm::vec3(0.0f);
					break;
				}
				//some step remains:
				remain *= (1.0f - time);
				//try to step over edge:
				glm::quat rotation;
				if (walkmesh->cross_edge(player.at, &end, &rotation)) {
					//stepped to a new triangle:
					player.at = end;
					//rotate step to follow surface:
					remain = rotation * remain;
				} else {
					//ran into a wall, bounce / slide along it:
					glm::vec3 const &a = walkmesh->vertices[player.at.indices.x];
					glm::vec3 const &b = walkmesh->vertices[player.at.indices.y];
					glm::vec3 const &c = walkmesh->vertices[player.at.indices.z];
					glm::vec3 along = glm::normalize(b-a);
					glm::vec3 normal = glm::normalize(glm::cross(b-a, c-a));
					glm::vec3 in = glm::cross(normal, along);

					//check how much 'remain' is pointing out of the triangle:
					float d = glm::dot(remain, in);
					if (d < 0.0f) {
						//bounce off of the wall:
						remain += (-1.25f * d) * in;
					} else {
						//if it's just pointing along the edge, bend slightly away from wall:
						remain += 0.01f * d * in;
					}
				}
			}

			if (remain != glm::vec3(0.0f)) {
				std::cout << "NOTE: code used full iteration budget for walking." << std::endl;
			}

			//update player's position to respect walking:
			player.transform->position = walkmesh->to_world_point(player.at);

			glm::vec3 const &a = walkmesh->vertices[player.at.indices.x];
			glm::vec3 const &b = walkmesh->vertices[player.at.indices.y];
			glm::vec3 const &c = walkmesh->vertices[player.at.indices.z];
			glm::vec3 normal = glm::cross(b-a, c-a);

			glm::vec3 player_pos = player.transform->position;
			glm::vec3 box_pos = box->position;
			glm::vec3 cone1_pos = cone1->position;
			glm::vec3 cone2_pos = cone2->position;
			glm::vec3 cylinder_long_pos = cylinder_long->position;
			glm::vec3 diff_box = box_pos - player_pos;
			glm::vec3 diff_cone1 = cone1_pos - player_pos;
			glm::vec3 diff_cone2 = cone2_pos - player_pos;
			glm::vec3 diff_cylinder_long = cylinder_long_pos - player_pos;

			if (glm::length(diff_box) < 0.5f ||
				(glm::length(diff_cone1) < 1.0f && glm::dot(normal, diff_cone1) >= 0.0f) ||
				(glm::length(diff_cone2) < 1.0f && glm::dot(normal, diff_cone2) >= 0.0f) ||
				(glm::length(diff_cylinder_long) < 1.5f && ((int) cylinder_long_total) % 180 >= 75 && ((int) cylinder_long_total) % 180 <= 105)) {
				player.at = old_at;
				player.transform->position = old_pos;
			}

			if (try_to_catch) {
				// glm::vec3 plain1 = ball1->position;
				// glm::vec3 plain2 = ball2->position;
				glm::vec3 ball1_pos = ball1->make_local_to_world() * glm::vec4(0.08715f, 0.3468f, -3.62069f, 1.0f);
				glm::vec3 ball2_pos = ball2->make_local_to_world() * glm::vec4(-0.19209f, 0.323112f, -3.19839f, 1.0f);
				glm::vec3 diff_ball1 = ball1_pos - player_pos;
				glm::vec3 diff_ball2 = ball2_pos - player_pos;
				if (glm::length(diff_ball1) < 3.0f) {
					ball1_caught = true;
					ball1->scale = glm::vec3(0.0f);
					Sound::play_3D(*sound_effect_load, 1.0f, ball1->position, 10.0f);
				}
				if (glm::length(diff_ball2) < 3.0f) {
					ball2_caught = true;
					ball2->scale = glm::vec3(0.0f);
					Sound::play_3D(*sound_effect_load, 1.0f, ball2->position, 10.0f);
				}
				try_to_catch = false;
			}

			{ //update player's rotation to respect local (smooth) up-vector:
				
				glm::quat adjust = glm::rotation(
					player.transform->rotation * glm::vec3(0.0f, 0.0f, 1.0f), //current up vector
					walkmesh->to_world_smooth_normal(player.at) //smoothed up vector at walk location
				);
				player.transform->rotation = glm::normalize(adjust * player.transform->rotation);
			}
		}

		{
			// other assets' movement
			float delta = elapsed * inf_speed;
			total_inf_translation += delta;
			if (total_inf_translation > 1.0f || total_inf_translation < -1.0f) {
				inf_speed = -inf_speed;
				total_inf_translation = 0.0f;
			}
			infinity->position = infinity->position + delta * glm::vec3(0.0f, 1.0f, 0.0f);

			float box_delta = elapsed * box_speed;
			up_and_down -= abs(box_delta);
			if (up_and_down < 0.0f) {
				box_speed = -box_speed;
				up_and_down = 1.5f;
			}
			box->position = box->position + box_delta * glm::vec3(0.0f, 1.0f, 0.0f);

			swing += elapsed / 10.0f;
			swing -= std::floor(swing);

			line1->rotation = line1_base_rotation * glm::angleAxis(
				glm::radians(5.0f * std::sin(swing * 2.0f * float(M_PI))),
				glm::vec3(1.0f, 0.0f, 0.0f)
			);
			line2->rotation = line2_base_rotation * glm::angleAxis(
				glm::radians(5.0f * std::sin(swing * 2.0f * float(M_PI))),
				glm::vec3(1.0f, 0.0f, 0.0f)
			);
			line3->rotation = line3_base_rotation * glm::angleAxis(
				glm::radians(5.0f * std::sin(swing * 2.0f * float(M_PI))),
				glm::vec3(1.0f, 0.0f, 0.0f)
			);

			ball1->rotation = ball1->rotation * glm::angleAxis(
				glm::radians(elapsed * 20.0f),
				glm::vec3(1.0f, 0.0f, 0.0f)
			);

			ball2->rotation = ball2->rotation * glm::angleAxis(
				glm::radians(elapsed * 50.0f),
				glm::vec3(0.0f, 1.0f, 0.0f)
			);

			float cone_delta = elapsed * cone_speed;
			cone_total -= abs(cone_delta);
			cone1->position = cone1->position + cone_delta * glm::vec3(0.0f, 1.0f, 0.0f);
			cone2->position = cone2->position - cone_delta * glm::vec3(0.0f, 1.0f, 0.0f);
			if (cone_total <= 0.0f) {
				cone_speed = -cone_speed;
				cone_total = 2.0f;
			}

			cylinder_long->rotation = cylinder_long->rotation * glm::angleAxis(
				glm::radians(elapsed * 50.0f),
				glm::vec3(1.0f, 0.0f, 0.0f)
			);
			cylinder_long_total += elapsed * 50.0f;
			if (cylinder_long_total > 360.0f) cylinder_long_total -= 360.0f;

		}

		{
			// text control
			if (!to_do) {
				update_window -= elapsed;
				if (update_window <= 0.0f) {
					cur_string_index += 1;
					to_do = true;
					update_window = 2.0f;
				}
			}
		}
	}

	//reset button press counters:
	left.downs = 0;
	right.downs = 0;
	up.downs = 0;
	down.downs = 0;
}

void PlayMode::draw(glm::uvec2 const &drawable_size) {
	//update camera aspect ratio for drawable:
	player.camera->aspect = float(drawable_size.x) / float(drawable_size.y);

	//set up light type and position for lit_color_texture_program:
	// TODO: consider using the Light(s) in the scene to do this
	glUseProgram(lit_color_texture_program->program);
	glUniform1i(lit_color_texture_program->LIGHT_TYPE_int, 1);
	glUniform3fv(lit_color_texture_program->LIGHT_DIRECTION_vec3, 1, glm::value_ptr(glm::vec3(0.0f, 0.0f,-1.0f)));
	glUniform3fv(lit_color_texture_program->LIGHT_ENERGY_vec3, 1, glm::value_ptr(glm::vec3(1.0f, 1.0f, 0.95f)));
	glUseProgram(0);

	glClearColor(0.5f, 0.5f, 0.8f, 1.0f);
	glClearDepth(1.0f); //1.0 is actually the default value to clear the depth buffer to, but FYI you can change it.
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS); //this is the default depth comparison function, but FYI you can change it.

	scene.draw(*player.camera);

	{ //use DrawLines to overlay some text:
		glDisable(GL_DEPTH_TEST);
		float aspect = float(drawable_size.x) / float(drawable_size.y);
		DrawLines lines(glm::mat4(
			1.0f / aspect, 0.0f, 0.0f, 0.0f,
			0.0f, 1.0f, 0.0f, 0.0f,
			0.0f, 0.0f, 1.0f, 0.0f,
			0.0f, 0.0f, 0.0f, 1.0f
		));

		constexpr float H = 0.09f;
		// lines.draw_text("Mouse motion looks; WASD moves; escape ungrabs mouse",
		// 	glm::vec3(-aspect + 0.1f * H, -1.0 + 0.1f * H, 0.0),
		// 	glm::vec3(H, 0.0f, 0.0f), glm::vec3(0.0f, H, 0.0f),
		// 	glm::u8vec4(0x00, 0x00, 0x00, 0x00));
		// float ofs = 2.0f / drawable_size.y;
		// lines.draw_text("Mouse motion looks; WASD moves; escape ungrabs mouse",
		// 	glm::vec3(-aspect + 0.1f * H + ofs, -1.0 + + 0.1f * H + ofs, 0.0),
		// 	glm::vec3(H, 0.0f, 0.0f), glm::vec3(0.0f, H, 0.0f),
		// 	glm::u8vec4(0xff, 0xff, 0xff, 0x00));
		if (to_do) {
			if (cur_string_index != 3) {
				lines.draw_text(hints[cur_string_index],
					glm::vec3(-aspect + 15.0f * H, -1.0 + 10.0f * H, 0.0),
					glm::vec3(H, 0.0f, 0.0f), glm::vec3(0.0f, H, 0.0f),
					glm::u8vec4(0xff, 0xff, 0xff, 0x00));
			} else {
				lines.draw_text(hints[cur_string_index],
					glm::vec3(-aspect + 0.1f * H, -1.0 + 0.1f * H, 0.0),
					glm::vec3(H, 0.0f, 0.0f), glm::vec3(0.0f, H, 0.0f),
					glm::u8vec4(0xff, 0xff, 0xff, 0x00));
			}
		}

	}
	GL_ERRORS();
}
