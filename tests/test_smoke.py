"""
Smoke tests for agentic workflow demo.
Tests core functionality with mocked external dependencies.
"""

import pytest
from unittest.mock import Mock, patch, MagicMock
import sys
import os

# Add src to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from main import parse_cli, config_logging, AgentState, create_workflow
from tools.compute import GlobusComputeWrapper
from sim_kernel import run_md_simulation


class TestCLI:
   """Test command line interface"""
   
   def test_parse_cli_defaults(self):
      """Test CLI parsing with default values"""
      with patch('sys.argv', ['main.py']):
         args = parse_cli()
         assert args.protein == "p53"
         assert args.log_level == "INFO"
         assert args.max_simulation_time == 60
   
   def test_parse_cli_custom_args(self):
      """Test CLI parsing with custom arguments"""
      test_args = [
         'main.py', 
         '--protein', 'insulin',
         '--model', 'gpt-4',
         '--log-level', 'DEBUG'
      ]
      with patch('sys.argv', test_args):
         args = parse_cli()
         assert args.protein == "insulin"
         assert args.model == "gpt-4"
         assert args.log_level == "DEBUG"


class TestSimulationKernel:
   """Test simulation kernel functionality"""
   
   def test_run_md_simulation_success(self):
      """Test successful simulation run"""
      params = {
         "protein": "p53",
         "timestep": 0.002,
         "temperature": 300,
         "steps": 1000
      }
      
      with patch('time.sleep'):  # Speed up test
         result = run_md_simulation(params)
      
      assert "error" not in result
      assert "final_energy" in result
      assert "rmsd" in result
      assert "stability_score" in result
      assert result["temperature"] == 300
      assert result["total_steps"] == 1000
   
   def test_run_md_simulation_reproducible(self):
      """Test that same protein gives reproducible results"""
      params = {
         "protein": "p53",
         "timestep": 0.002,
         "temperature": 300,
         "steps": 1000
      }
      
      with patch('time.sleep'):
         result1 = run_md_simulation(params)
         result2 = run_md_simulation(params)
      
      assert result1["final_energy"] == result2["final_energy"]
      assert result1["rmsd"] == result2["rmsd"]


class TestGlobusComputeWrapper:
   """Test Globus Compute wrapper"""
   
   def test_init_with_endpoint_id(self):
      """Test initialization with explicit endpoint ID"""
      # Use environment variable or fallback UUID for testing
      test_uuid = os.getenv("GC_ENDPOINT_ID", "12b14478-2b82-4d1a-8588-cc3a956735d1")
      with patch('globus_compute_sdk.Client'), \
           patch('globus_compute_sdk.Executor'):
         wrapper = GlobusComputeWrapper(test_uuid)
         assert wrapper.endpoint_id == test_uuid
   
   def test_init_no_endpoint_raises(self):
      """Test that missing endpoint ID raises ValueError"""
      with patch.dict(os.environ, {}, clear=True):
         with pytest.raises(ValueError, match="endpoint ID not provided"):
            GlobusComputeWrapper()
   
   @patch('globus_compute_sdk.Client')
   @patch('globus_compute_sdk.Executor')
   def test_submit_simulation_success(self, mock_executor_class, mock_client_class):
      """Test successful simulation submission"""
      # Mock the future object
      mock_future = Mock()
      mock_future.task_id = "task-123"
      mock_future.done.return_value = True
      mock_future.result.return_value = {
         "final_energy": -12500.0,
         "rmsd": 1.5,
         "stability_score": 0.8
      }
      
      # Mock the executor instance
      mock_executor = Mock()
      mock_executor.submit.return_value = mock_future
      mock_executor_class.return_value = mock_executor
      
      # Use environment variable or fallback UUID for testing
      test_uuid = os.getenv("GC_ENDPOINT_ID", "12b14478-2b82-4d1a-8588-cc3a956735d1")
      wrapper = GlobusComputeWrapper(test_uuid)
      
      # Override the executor instance with our mock
      wrapper.executor = mock_executor
      
      params = {"protein": "p53", "steps": 1000}
      result = wrapper.submit_simulation(params)
      
      assert result["status"] == "completed"
      assert result["task_id"] == "task-123"
      assert "final_energy" in result
   
   @patch('globus_compute_sdk.Client')
   @patch('globus_compute_sdk.Executor')
   def test_submit_simulation_timeout(self, mock_executor_class, mock_client_class):
      """Test simulation timeout handling"""
      # Mock future that never completes
      mock_future = Mock()
      mock_future.task_id = "task-456"
      mock_future.done.return_value = False
      
      mock_executor = Mock()
      mock_executor.submit.return_value = mock_future
      mock_executor_class.return_value = mock_executor
      
      # Use environment variable or fallback UUID for testing
      test_uuid = os.getenv("GC_ENDPOINT_ID", "12b14478-2b82-4d1a-8588-cc3a956735d1")
      wrapper = GlobusComputeWrapper(test_uuid)
      
      # Override the executor instance with our mock
      wrapper.executor = mock_executor
      
      params = {"protein": "p53", "steps": 1000}
      result = wrapper.submit_simulation(params, timeout=1)  # Very short timeout
      
      assert result["status"] == "timeout"
      assert "timed out" in result["error"]


class TestWorkflow:
   """Test LangGraph workflow"""
   
   @patch('main.ChatOpenAI')
   @patch('main.GlobusComputeWrapper')
   def test_workflow_creation(self, mock_gc, mock_llm):
      """Test that workflow can be created"""
      workflow = create_workflow()
      assert workflow is not None
   
   @patch('main.ChatOpenAI')
   @patch('main.GlobusComputeWrapper')
   def test_llm_analysis_node(self, mock_gc, mock_llm):
      """Test LLM analysis node"""
      from main import llm_analysis_node
      from langchain_core.messages import AIMessage
      
      # Mock LLM response as a proper LangChain message
      mock_response = AIMessage(content="Mock analysis response for p53 protein")
      mock_llm_instance = Mock()
      mock_llm_instance.invoke.return_value = mock_response
      mock_llm.return_value = mock_llm_instance
      
      state = AgentState(
         messages=[],
         protein="p53",
         analysis_request="",
         simulation_params={},
         simulation_result={},
         final_report=""
      )
      
      result_state = llm_analysis_node(state)
      
      assert "simulation_params" in result_state
      assert result_state["simulation_params"]["protein"] == "p53"
      assert len(result_state["messages"]) > 0
   
   @patch('main.GlobusComputeWrapper')
   def test_simulation_node(self, mock_gc_class):
      """Test simulation node"""
      from main import simulation_node
      
      # Mock GC wrapper
      mock_gc_instance = Mock()
      mock_gc_instance.submit_simulation.return_value = {
         "status": "completed",
         "final_energy": -12000,
         "rmsd": 1.2
      }
      mock_gc_class.return_value = mock_gc_instance
      
      state = AgentState(
         messages=[],
         protein="insulin",
         analysis_request="",
         simulation_params={"protein": "insulin", "steps": 5000},
         simulation_result={},
         final_report=""
      )
      
      result_state = simulation_node(state)
      
      assert result_state["simulation_result"]["status"] == "completed"
      assert "final_energy" in result_state["simulation_result"]
   
   def test_report_node_success(self):
      """Test report generation for successful simulation"""
      from main import report_node
      
      state = AgentState(
         messages=[],
         protein="p53",
         analysis_request="",
         simulation_params={"temperature": 300, "timestep": 0.002, "steps": 10000},
         simulation_result={
            "final_energy": -12500.0,
            "rmsd": 1.3,
            "stability_score": 0.85
         },
         final_report=""
      )
      
      result_state = report_node(state)
      
      assert "COMPLETED" in result_state["final_report"]
      assert "p53" in result_state["final_report"]
      assert "-12500.0" in result_state["final_report"]
   
   def test_report_node_failure(self):
      """Test report generation for failed simulation"""
      from main import report_node
      
      state = AgentState(
         messages=[],
         protein="p53",
         analysis_request="",
         simulation_params={},
         simulation_result={"error": "Endpoint unreachable"},
         final_report=""
      )
      
      result_state = report_node(state)
      
      assert "FAILED" in result_state["final_report"]
      assert "Endpoint unreachable" in result_state["final_report"]


if __name__ == "__main__":
   pytest.main([__file__, "-v"]) 