import mujoco
import numpy as np

from .lcm2mujuco_bridge import Lcm2MujocoBridge
from arc_bridge.lcm_msgs import bdx_droid_state_t, bdx_droid_control_t
from arc_bridge.utils import *


class BdxDroidBridge(Lcm2MujocoBridge):
    def __init__(self, mj_model, mj_data, config):
        super().__init__(mj_model, mj_data, config)

        self.right_foot_name = "right_foot"
        self.left_foot_name = "left_foot"
        self.height_init = 1.0

    def parse_robot_specific_low_state(self):
        # """Add robot-specific state information to low_state message"""
        # Example: Add inertia matrix and bias forces
        temp_inertia_mat = np.zeros((self.mj_model.nv, self.mj_model.nv))
        mujoco.mj_fullM(self.mj_model, temp_inertia_mat, self.mj_data.qM)
        self.low_state.inertia_mat = temp_inertia_mat.tolist()
        self.low_state.bias_force = self.mj_data.qfrc_bias.tolist()

        # tau = Mddq + C + G, C+G is the bias forces
        self.update_kinematics()
    
    def update_kinematics(self):
        dq = np.zeros((self.mj_model.nv, )) 

        # Right foot Jacobian
        # translational J of right-foot site, 3xnv
        right_foot_id = mujoco.mj_name2id(self.mj_model, mujoco.mjtObj.mjOBJ_SITE, self.right_foot_name)

        # J_linear_foot_R (aka linear velocity)
        right_foot_pos = self.mj_data.site_xpos[right_foot_id] # current world position
        J_foot_R = np.zeros((3, self.mj_model.nv))
        mujoco.mj_jacSite(self.mj_model, self.mj_data, J_foot_R, None, right_foot_id)
        
        # dJ/dq R
        dJ_foot_R = np.zeros((3, self.mj_model.nv))
        mujoco.mj_jacDot(self.mj_model, self.mj_data, dJ_foot_R, None, right_foot_pos, right_foot_id) # needs body site? not right_foot_id
        dJdq_foot_R = dJ_foot_R @ dq 

        # Left foot Jacobian. why getting it at site instead of body? 
        # translational J of left-foot site, 3xnv
        left_foot_id = mujoco.mj_name2id(self.mj_model, mujoco.mjtObj.mjOBJ_SITE, self.left_foot_name) 
        
        # J_linear_foot_L (aka linear velocity)
        left_foot_pos = self.mj_data.site_xpos[left_foot_id]
        J_foot_L = np.zeros((3, self.mj_model.nv))
        mujoco.mj_jacSite(self.mj_model, self.mj_data, J_foot_L, None, left_foot_id)

        # dJ/dq L
        dJ_foot_L = np.zeros((3, self.mj_model.nv))
        mujoco.mj_jacDot(self.mj_model, self.mj_data, dJ_foot_L, None, left_foot_pos, left_foot_id) # needs body site? not right_foot_id
        dJdq_foot_L = dJ_foot_L @ dq 

        # RIGHT NOW BASED ON BIPEDPOINTFOOT LCM TYPE -> change to rabbit? 
        # wouldn't BIPEDPOINTFOOT or BIPEDLINEFOOT would better for us?
        # should pf be vstack?
        # gc - ground contact
        p_gc = np.hstack((right_foot_pos, left_foot_pos)) # concatenate positions, 6x1
        J_gc = np.vstack((J_foot_R, J_foot_L)) # stack Jacobians, 6xnv
        # dJdq_gc = np.vstack((dJdq_foot_R, dJdq_foot_L)) # 6x1 ?
        dJdq_gc = np.hstack((dJdq_foot_R,dJdq_foot_L))
        # send through LCM
        self.low_state.p_gc = p_gc.tolist()
        self.low_state.J_gc = J_gc.tolist()
        self.low_state.dJdq_gc = dJdq_gc.tolist()

    def lcm_state_handler(self, channel, data):
        if self.mj_data == None:
            return
        msg = bdx_droid_state_t.decode(data)
